"""
EVE - Network Interceptor & Cryptanalysis Station
Captures raw TCP packets between Alice and Bob, displays hex dumps, runs ECDLP attacks to recover keys.
"""
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import eventlet
eventlet.monkey_patch()

import sys
import os
import time
import subprocess
import random
import json
import threading
from pathlib import Path
from flask import Flask, render_template, request
from flask_socketio import SocketIO, emit
from scapy.all import sniff, TCP, IP, Raw

# Project Root Setup
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(PROJECT_ROOT)

try:
    from utils import load_input, EllipticCurve
    from web_interface.crypto_manager import CryptoManager
except ImportError as e:
    print(f"IMPORT ERROR: {e}")

app = Flask(__name__)
socketio = SocketIO(app, cors_allowed_origins="*", async_mode='eventlet')

# --- GLOBAL STATE ---
interceptor_state = {
    "running": False,
    "packets": [],
    "alice_pub": None,
    "bob_pub": None,
    "curve_params": None,
    "alice_priv": None,  # Cracked Alice's private key
    "bob_priv": None,    # Cracked Bob's private key
    "aes_key": None,
    "attack_process": None,
    "attack_running": False
}

@app.route('/')
def home():
    return render_template('eve_interceptor.html')

@app.route('/packet_sent', methods=['POST'])
def receive_packet():
    """Receive packet info directly from Alice/Bob (Bypasses Scapy need for root)"""
    import re
    data = request.get_json()
    
    if data:
        # Store public keys if this is a handshake
        decoded_text = data.get('decoded', '')
        print(f"[EVE DEBUG] Received from {data.get('src_label')}: {decoded_text}")
        
        if 'PubKey' in decoded_text or 'pubkey' in decoded_text.lower():
            if data['src_label'] == 'ALICE':
                # Parse public key from decoded string
                match = re.search(r'\((\d+), (\d+)\)', decoded_text)
                if match:
                    x, y = int(match.group(1)), int(match.group(2))
                    interceptor_state["alice_pub"] = (x, y)
                    print(f"[EVE] ✓ Captured Alice's public key: ({x}, {y})")
                else:
                    print(f"[EVE] ✗ Failed to parse Alice's coordinates")
                
                # Extract bit length and case from format: "PubKey[35bit,case1]: (x, y)"
                bit_case_match = re.search(r'\[(\d+)bit,case(\d+)\]', decoded_text)
                if bit_case_match:
                    interceptor_state["alice_bits"] = int(bit_case_match.group(1))
                    interceptor_state["alice_case"] = int(bit_case_match.group(2))
                    print(f"[EVE] ✓ Detected Alice's curve: {interceptor_state['alice_bits']}-bit, case {interceptor_state['alice_case']}")
                else:
                    print(f"[EVE] ✗ Failed to parse Alice's curve metadata")
                    
            elif data['src_label'] == 'BOB':
                match = re.search(r'\((\d+), (\d+)\)', decoded_text)
                if match:
                    x, y = int(match.group(1)), int(match.group(2))
                    interceptor_state["bob_pub"] = (x, y)
                    print(f"[EVE] ✓ Captured Bob's public key: ({x}, {y})")
                else:
                    print(f"[EVE] ✗ Failed to parse Bob's coordinates")
                
                bit_case_match = re.search(r'\[(\d+)bit,case(\d+)\]', decoded_text)
                if bit_case_match:
                    interceptor_state["bob_bits"] = int(bit_case_match.group(1))
                    interceptor_state["bob_case"] = int(bit_case_match.group(2))
                    print(f"[EVE] ✓ Detected Bob's curve: {interceptor_state['bob_bits']}-bit, case {interceptor_state['bob_case']}")
                else:
                    print(f"[EVE] ✗ Failed to parse Bob's curve metadata")
        
        # Broadcast to all connected clients
        interceptor_state["packets"].append(data)
        socketio.emit('packet_captured', data)
        
        # Debug: Show captured keys
        if interceptor_state.get("alice_pub"):
            print(f"[EVE] Status: Alice key captured ✓")
        if interceptor_state.get("bob_pub"):
            print(f"[EVE] Status: Bob key captured ✓")
    
    return '', 204

# --- PACKET SNIFFER (Optional fallback) ---
def packet_handler(pkt):
    if not interceptor_state["running"]: return
    if pkt.haslayer(TCP) and pkt.haslayer(Raw):
        # (Scapy logic kept for completeness, but /packet_sent is preferred for demo stability)
        pass 

@socketio.on('start_capture')
def handle_start():
    interceptor_state["running"] = True
    emit('log', {'text': "[EVE] Capture ACTIVE. Listening for Alice/Bob traffic..."})

@socketio.on('stop_capture')
def handle_stop():
    interceptor_state["running"] = False
    # Kill any running attack
    if interceptor_state["attack_running"] and interceptor_state["attack_process"]:
        try:
            interceptor_state["attack_process"].terminate()
        except: pass
        interceptor_state["attack_running"] = False
    
    interceptor_state["alice_priv"] = None
    interceptor_state["bob_priv"] = None
    interceptor_state["aes_key"] = None
    emit('log', {'text': "[EVE] Capture STOPPED."})
    emit('attack_reset')

@socketio.on('run_attack')
def handle_attack(data):
    algo = data.get('algo', 'PollardRho')
    leak = int(data.get('leak', 0))
    mode = data.get('mode', 'basic')
    interval_width = data.get('interval_width', None)
    gaussian = data.get('gaussian', False)
    packet_index = data.get('packet_index', None)
    
    # Use Alice's actual curve parameters from handshake
    if not interceptor_state.get("alice_bits") or not interceptor_state.get("alice_case"):
        emit('log', {'text': "[ERROR] No handshake captured yet! Wait for Alice/Bob to connect."})
        return
    
    curve_bits = interceptor_state["alice_bits"]
    case_num = interceptor_state["alice_case"]
    
    emit('log', {'text': f"[INFO] Using Alice's actual curve: {curve_bits}-bit, case {case_num}"})
    
    # 1. Determine Target Public Key
    target_pub_key = None
    
    # Try to get from specific packet
    if packet_index is not None and 0 <= packet_index < len(interceptor_state["packets"]):
        packet = interceptor_state["packets"][packet_index]
        if packet.get('decoded', '').startswith('PubKey'):
            import re
            match = re.search(r'\((\d+), (\d+)\)', packet['decoded'])
            if match:
                target_pub_key = (int(match.group(1)), int(match.group(2)))
                emit('log', {'text': f"[ATTACK] Targeted Packet #{packet_index + 1}"})

    # Fallback to captured Alice key
    if not target_pub_key:
        if not interceptor_state["alice_pub"]:
            emit('log', {'text': "[ERROR] No Public Key found! Wait for handshake."})
            return
        target_pub_key = interceptor_state["alice_pub"]
        emit('log', {'text': "[ATTACK] Targeting Alice's Public Key"})

    interceptor_state["attack_running"] = True
    
    # 2. Prepare Test Case Files
    case_path = Path(PROJECT_ROOT) / 'test_cases' / f'{curve_bits}bit' / f'case_{case_num}.txt'
    if not case_path.exists():
        emit('log', {'text': f"[ERROR] Case file not found: {case_path}"})
        interceptor_state["attack_running"] = False
        return

    try:
        p, a, b, G, n, _ = load_input(case_path) #type:ignore
        interceptor_state["curve_params"] = {'p': p, 'a': a, 'b': b, 'G': G, 'n': n}
        
        # Write Temp Case for the script to read
        temp_case = Path(PROJECT_ROOT) / 'web_interface' / 'eve_temp.txt'
        with open(temp_case, 'w') as f:
            f.write(f"{p}\n{a} {b}\n{G[0]} {G[1]}\n{n}\n{target_pub_key[0]} {target_pub_key[1]}\n") #type:ignore

        # 3. Prepare Command
        scripts = {
            'BruteForce': 'BruteForce/main_optimized.py',
            'BruteForce_Bonus': 'BruteForce/bonus.py',
            'BabyStep': 'BabyStep/main_optimized.py',
            'PollardRho': 'PollardRho/main_optimized.py',
            'PohligHellman': 'PohligHellman/main_optimized.py',
            'LasVegas': 'LasVegas/main_optimized.py',
            'LasVegas_Gaussian': 'LasVegas/bonus.py',
        }
        
        # Choose script based on mode
        if gaussian and algo == 'LasVegas':
            script_key = 'LasVegas_Gaussian'
            emit('log', {'text': "[MODE] Using Las Vegas Gaussian Estimation"})
        elif (leak > 0 or interval_width) and algo == 'BruteForce':
            script_key = 'BruteForce_Bonus'
            emit('log', {'text': f"[MODE] Using BruteForce with {'LSB leak' if leak > 0 else 'bounded interval'}"})
        else:
            script_key = algo
        
        script_rel = scripts.get(script_key, 'PollardRho/main_optimized.py')
        script_path = Path(PROJECT_ROOT) / script_rel
        
        if not script_path.exists():
            emit('log', {'text': f"[ERROR] Script not found: {script_rel}"})
            return

        venv_python = Path(PROJECT_ROOT) / 'venv' / 'bin' / 'python3'
        python_exe = str(venv_python) if venv_python.exists() else sys.executable
        
        cmd = [python_exe, "-u", str(script_path), str(temp_case)]
        
        # Add mode-specific arguments
        if leak > 0 and script_key == 'BruteForce_Bonus':
            cmd.extend(["--leak-bits", str(leak)])
            emit('log', {'text': f"[ATTACK] Using {leak}-bit LSB leak"})
        elif interval_width and script_key == 'BruteForce_Bonus':
            cmd.extend(["--interval-width", str(interval_width)])
            emit('log', {'text': f"[ATTACK] Searching bounded interval: ±{interval_width//2}"})
        elif leak > 0 and script_key != 'BruteForce_Bonus':
            # Standard algorithms with leak support
            cmd.extend(["--leak-bits", str(leak)])

        emit('log', {'text': f"[EXEC] Running {algo} on {curve_bits}-bit curve..."})

        # 4. Run Process with Real-Time Output
        def execution_thread():
            solution_found = False
            private_key = None
            
            try:
                # IMPORTANT: Unbuffered output environment
                env = os.environ.copy()
                env["PYTHONUNBUFFERED"] = "1"
                
                proc = subprocess.Popen(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, cwd=PROJECT_ROOT, env=env, bufsize=0
                )
                interceptor_state["attack_process"] = proc
                
                # Character-by-character reading loop
                buffer = []
                while True:
                    if not interceptor_state["attack_running"]:
                        proc.terminate()
                        socketio.emit('log', {'text': "[INFO] Attack aborted."})
                        break
                    
                    char = proc.stdout.read(1) #type:ignore
                    if not char and proc.poll() is not None:
                        break
                    if not char: continue
                    
                    # Flush on newline OR carriage return (progress bars)
                    if char in ['\r', '\n']:
                        line = "".join(buffer).strip()
                        if line:
                            socketio.emit('log', {'text': f"  {line}"})
                            # Check for success
                            if "Solution: d =" in line:
                                solution_found = True
                                try: private_key = int(line.split("d =")[1].strip().split()[0])
                                except: pass
                        buffer = []
                    else:
                        buffer.append(char)
                
                # Flush remaining
                if buffer: 
                    line = "".join(buffer).strip()
                    if line: socketio.emit('log', {'text': f"  {line}"})
                
                proc.wait()
                
                if solution_found and private_key:
                    # Determine whose key we cracked
                    if target_pub_key == interceptor_state.get("alice_pub"):
                        interceptor_state["alice_priv"] = private_key
                        socketio.emit('log', {'text': f"[SUCCESS] ✓ CRACKED Alice's Private Key: d_alice = {private_key}"})
                        socketio.emit('crack_success', {'alice_key': str(private_key), 'target': 'ALICE'})
                    elif target_pub_key == interceptor_state.get("bob_pub"):
                        interceptor_state["bob_priv"] = private_key
                        socketio.emit('log', {'text': f"[SUCCESS] ✓ CRACKED Bob's Private Key: d_bob = {private_key}"})
                        socketio.emit('crack_success', {'bob_key': str(private_key), 'target': 'BOB'})
                    socketio.emit('log', {'text': f"[INFO] You can now decrypt {('Alice→Bob' if target_pub_key == interceptor_state.get('alice_pub') else 'Bob→Alice')} messages."})
                else:
                    if interceptor_state["attack_running"]:
                        socketio.emit('log', {'text': "[FAILED] Attack finished without finding key."})

            except Exception as e:
                socketio.emit('log', {'text': f"[ERROR] Runtime exception: {e}"})
            finally:
                interceptor_state["attack_running"] = False
        
        socketio.start_background_task(execution_thread)

    except Exception as e:
        emit('log', {'text': f"[ERROR] {e}"})
        interceptor_state["attack_running"] = False

@socketio.on('decrypt_packet')
def handle_decrypt(data):
    """Decrypt a single packet using provided private key d"""
    idx = data.get('packet_index')
    private_key_d = data.get('private_key')
    
    if not private_key_d:
        emit('log', {'text': "[ERROR] No private key provided! Enter the cracked private key."})
        return
    
    try:
        private_key_d = int(private_key_d)
    except:
        emit('log', {'text': "[ERROR] Invalid private key format. Must be a number."})
        return
    
    if not interceptor_state.get("alice_pub") or not interceptor_state.get("bob_pub"):
        emit('log', {'text': "[ERROR] Public keys not captured. Wait for handshake."})
        return
        
    pkt = interceptor_state["packets"][idx]
    sender = pkt.get('src_label', 'UNKNOWN')
    
    try:
        emit('log', {'text': f"═══════════════════════════════════════════════════"})
        emit('log', {'text': f"[DECRYPT] Packet from {sender}"})
        
        # Determine which key to use based on sender
        params = interceptor_state["curve_params"]
        curve = EllipticCurve(params['a'], params['b'], params['p']) #type:ignore
        
        if sender == 'ALICE':
            # Alice sent this: decrypt using d_alice × Q_bob
            emit('log', {'text': f"[INFO] Using d_alice = {private_key_d} × Q_bob"})
            receiver_pub = interceptor_state["bob_pub"]
        elif sender == 'BOB':
            # Bob sent this: decrypt using d_bob × Q_alice
            emit('log', {'text': f"[INFO] Using d_bob = {private_key_d} × Q_alice"})
            receiver_pub = interceptor_state["alice_pub"]
        else:
            emit('log', {'text': f"[ERROR] Unknown sender: {sender}"})
            return
        
        emit('log', {'text': f"[INFO] Sender's private key (cracked): d = {private_key_d}"})
        emit('log', {'text': f"[INFO] Receiver's public key: Q = {receiver_pub}"})
        
        # Derive shared secret: sender_priv × receiver_pub
        shared = curve.scalar_multiply(private_key_d, receiver_pub)
        
        if not shared:
            emit('log', {'text': "[ERROR] Shared secret is Infinity!"})
            return
        
        emit('log', {'text': f"[INFO] Shared Secret: S = ({shared[0]}, {shared[1]})"})
        
        # Derive AES key
        aes_key = CryptoManager.derive_key(shared[0]) #type:ignore
        emit('log', {'text': f"[INFO] AES Key: {aes_key.hex()}"})
        
        # Parse encrypted packet
        hex_data = pkt['hex']
        data_bytes = bytes.fromhex(hex_data)
        
        nl = int.from_bytes(data_bytes[0:4], 'big')
        nonce = data_bytes[4:4+nl]
        cl = int.from_bytes(data_bytes[4+nl:8+nl], 'big')
        ct = data_bytes[8+nl:8+nl+cl]
        
        emit('log', {'text': f"[ENCRYPTED] Nonce: {nonce.hex()}"})
        emit('log', {'text': f"[ENCRYPTED] Ciphertext: {ct.hex()[:64]}... ({len(ct)} bytes)"})
        
        # Decrypt
        plain = CryptoManager.decrypt(aes_key, nonce.hex(), ct.hex()) #type:ignore
        
        emit('log', {'text': f"[DECRYPTED] ✓ {sender}: \"{plain}\""})
        emit('log', {'text': f"═══════════════════════════════════════════════════"})
        
        socketio.emit('decrypted_message', {
            'packet_index': idx,
            'plaintext': plain,
            'sender': sender,
            'encrypted_hex': ct.hex(),
            'nonce': nonce.hex()
        })
    except Exception as e:
        import traceback
        emit('log', {'text': f"[ERROR] Decryption failed: {e}"})
        emit('log', {'text': f"[ERROR] Traceback: {traceback.format_exc()[:200]}"})

@socketio.on('decrypt_all')
def decrypt_all():
    """Decrypt everything in buffer"""
    if not interceptor_state.get("aes_key"): return
    count = 0
    for i, pkt in enumerate(interceptor_state["packets"]):
        if 'AES' in pkt.get('decoded', ''):
            handle_decrypt({'packet_index': i})
            count += 1
    emit('log', {'text': f"[INFO] Decrypted {count} messages."})

if __name__ == '__main__':
    print("[EVE] Server running on http://127.0.0.1:5500")
    socketio.run(app, debug=True, port=5500)