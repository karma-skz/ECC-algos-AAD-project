"""ALICE - Sender"""
import socket
import threading
import time
import random
import sys
import argparse
import requests
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from utils import load_input, EllipticCurve
from web_interface.crypto_manager import CryptoManager

ALICE_PORT = 5001
EVE_URL = "http://localhost:5500/packet_sent"

class AliceServer:
    def __init__(self, test_case_path):
        self.test_case_path = Path(test_case_path)
        if not self.test_case_path.exists():
            raise FileNotFoundError(f"Test case not found: {test_case_path}")
        
        self.sock = None
        self.conn = None
        self.aes_key = None
        self.setup_ecc()

    def setup_ecc(self):
        p, a, b, G, n, Q_test = load_input(self.test_case_path)
        self.curve = EllipticCurve(a, b, p)
        
        # Extract bit length and case from path: test_cases/35bit/case_1.txt
        import re
        path_str = str(self.test_case_path)
        bit_match = re.search(r'(\d+)bit', path_str)
        case_match = re.search(r'case_(\d+)', path_str)
        self.bits = int(bit_match.group(1)) if bit_match else 30
        self.case = int(case_match.group(1)) if case_match else 1
        
        # Use the test case's private key from answer file
        answer_file = self.test_case_path.parent / f'answer_{self.case}.txt'
        if answer_file.exists():
            with open(answer_file) as af:
                self.priv = int(af.read().strip())
        else:
            # Fallback: use a small known value that's crackable
            self.priv = random.randint(1, min(10000, n-1))
        
        self.pub = self.curve.scalar_multiply(self.priv, G)
        
        print(f"[ALICE] {self.bits}-bit ECC (Case {self.case})")
        print(f"[ALICE] Private key: d = {self.priv}")
        print(f"[ALICE] Public key: Q = ({self.pub[0]}, {self.pub[1]})")#type:ignore
        print(f"[ALICE] Order: n = {n}")

    def notify_eve(self, label, data_hex, decoded):
        try:
            requests.post(EVE_URL, json={
                "time": time.strftime("%H:%M:%S"),
                "src": f"127.0.0.1:{ALICE_PORT}", "dst": "127.0.0.1:5002",
                "src_label": "ALICE", "dst_label": "BOB",
                "proto": label, "len": len(data_hex)//2, "hex": data_hex, "decoded": decoded
            }, timeout=0.1)
        except: pass

    def run(self):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.bind(("0.0.0.0", ALICE_PORT))
        self.sock.listen(1)
        print(f"[ALICE] Listening on {ALICE_PORT}...")
        
        try:
            self.conn, addr = self.sock.accept()
            print(f"[ALICE] Connected to {addr}")
            
            # 1. Send PubKey
            pk_str = f"{self.pub[0]},{self.pub[1]}".encode()#type:ignore
            self.conn.sendall(len(pk_str).to_bytes(4, 'big') + pk_str)
            self.notify_eve("ECDH", pk_str.hex(), f"PubKey[{self.bits}bit,case{self.case}]: ({self.pub[0]}, {self.pub[1]})")#type:ignore
            
            
            # 2. Receive PubKey
            l = int.from_bytes(self.conn.recv(4), 'big')
            bob_pk_str = self.conn.recv(l).decode()
            bx, by = map(int, bob_pk_str.split(','))
            
            # 3. Derive Secret
            S = self.curve.scalar_multiply(self.priv, (bx, by))
            self.aes_key = CryptoManager.derive_key(S[0])#type:ignore
            print(f"[ALICE] Secure Channel Established. AES Key Derived.")
            print(f"[ALICE DEBUG] Shared Secret: S = ({S[0]}, {S[1]})")#type:ignore
            print(f"[ALICE DEBUG] AES Key: {self.aes_key.hex()}")
            
            # 4. Chat Loop
            threading.Thread(target=self.rx_loop, daemon=True).start()
            while True:
                msg = input("Alice> ")
                if not msg: continue
                
                enc = CryptoManager.encrypt(self.aes_key, msg)
                nb = bytes.fromhex(enc['nonce'])
                cb = bytes.fromhex(enc['ciphertext'])
                
                # Packet: [nonce_len 4][nonce][ct_len 4][ct]
                pkt = len(nb).to_bytes(4,'big') + nb + len(cb).to_bytes(4,'big') + cb
                self.conn.sendall(pkt)
                self.notify_eve("TCP/TLS", pkt.hex(), "AES-256 Encrypted")
        except KeyboardInterrupt:
            print("\n[ALICE] Shutting down...")
        except Exception as e:
            print(f"\n[ALICE] Error: {e}")
        finally:
            if self.conn:
                try:
                    self.conn.close()
                except:
                    pass
            if self.sock:
                try:
                    self.sock.close()
                except:
                    pass
            print("[ALICE] Disconnected.")

    def rx_loop(self):
        try:
            while True:
                nl = int.from_bytes(self.conn.recv(4), 'big')#type:ignore
                if nl == 0: break
                n = self.conn.recv(nl).hex()#type:ignore
                cl = int.from_bytes(self.conn.recv(4), 'big')#type:ignore
                c = self.conn.recv(cl).hex()#type:ignore
                
                txt = CryptoManager.decrypt(self.aes_key, n, c) #type:ignore
                print(f"\r[Bob]: {txt}\nAlice> ", end="")
        except Exception as e:
            print(f"\n[ALICE RX] Connection closed: {e}")
            pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Alice Server - ECDH Sender')
    parser.add_argument('testcase', type=str, help='Path to test case file (e.g., test_cases/35bit/case_1.txt)')
    args = parser.parse_args()
    
    AliceServer(args.testcase).run()