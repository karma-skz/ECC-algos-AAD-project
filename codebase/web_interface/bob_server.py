"""BOB - Receiver"""
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

BOB_PORT = 5002
EVE_URL = "http://localhost:5500/packet_sent"

class BobServer:
    def __init__(self, test_case_path):
        self.test_case_path = Path(test_case_path)
        if not self.test_case_path.exists():
            raise FileNotFoundError(f"Test case not found: {test_case_path}")
        self.setup_ecc()

    def setup_ecc(self):
        p, a, b, G, n, Q_test = load_input(self.test_case_path)
        self.curve = EllipticCurve(a, b, p)
        
        # Extract bit length and case from path
        import re
        path_str = str(self.test_case_path)
        bit_match = re.search(r'(\d+)bit', path_str)
        case_match = re.search(r'case_(\d+)', path_str)
        self.bits = int(bit_match.group(1)) if bit_match else 30
        self.case = int(case_match.group(1)) if case_match else 1
        
        # Use a different private key than Alice (but still crackable)
        self.priv = random.randint(1, min(10000, n-1))
        self.pub = self.curve.scalar_multiply(self.priv, G)
        
        print(f"[BOB] {self.bits}-bit ECC (Case {self.case})")
        print(f"[BOB] Private key: d = {self.priv}")
        print(f"[BOB] Public key: Q = ({self.pub[0]}, {self.pub[1]})")#type:ignore

    def notify_eve(self, label, data_hex, decoded):
        try:
            requests.post(EVE_URL, json={
                "time": time.strftime("%H:%M:%S"),
                "src": f"127.0.0.1:{BOB_PORT}", "dst": "127.0.0.1:5001",
                "src_label": "BOB", "dst_label": "ALICE",
                "proto": label, "len": len(data_hex)//2, "hex": data_hex, "decoded": decoded
            }, timeout=0.1)
        except: pass

    def run(self):
        print("[BOB] Connecting to Alice...")
        while True:
            try:
                self.conn = socket.create_connection(("127.0.0.1", 5001), timeout=5)
                break
            except:
                print("[BOB] Waiting for Alice...")
                time.sleep(1)
        
        try:
            # 1. Receive Alice Pub
            l = int.from_bytes(self.conn.recv(4), 'big')
            alice_pk_str = self.conn.recv(l).decode()
            ax, ay = map(int, alice_pk_str.split(','))
            
            # 2. Send Bob Pub
            pk_str = f"{self.pub[0]},{self.pub[1]}".encode()#type:ignore
            self.conn.sendall(len(pk_str).to_bytes(4, 'big') + pk_str)
            self.notify_eve("ECDH", pk_str.hex(), f"PubKey[{self.bits}bit,case{self.case}]: ({self.pub[0]}, {self.pub[1]})")#type:ignore
            
            # 3. Derive
            S = self.curve.scalar_multiply(self.priv, (ax, ay))
            self.aes_key = CryptoManager.derive_key(S[0])#type:ignore
            print(f"[BOB] Secure Channel Established.")
            print(f"[BOB DEBUG] Shared Secret: S = ({S[0]}, {S[1]})")#type:ignore
            print(f"[BOB DEBUG] AES Key: {self.aes_key.hex()}")#type:ignore
            
            threading.Thread(target=self.rx_loop, daemon=True).start()
            while True:
                msg = input("Bob> ")
                if not msg: continue
                enc = CryptoManager.encrypt(self.aes_key, msg)
                nb = bytes.fromhex(enc['nonce'])
                cb = bytes.fromhex(enc['ciphertext'])
                
                # Packet: [nonce_len 4][nonce][ct_len 4][ct]
                pkt = len(nb).to_bytes(4,'big') + nb + len(cb).to_bytes(4,'big') + cb
                self.conn.sendall(pkt)
                self.notify_eve("TCP/TLS", pkt.hex(), "AES-256 Encrypted")
        except KeyboardInterrupt:
            print("\n[BOB] Shutting down...")
        except Exception as e:
            print(f"\n[BOB] Error: {e}")
        finally:
            if self.conn:
                try:
                    self.conn.close()
                except:
                    pass
            print("[BOB] Disconnected.")

    def rx_loop(self):
        try:
            while True:
                nl = int.from_bytes(self.conn.recv(4), 'big')
                if nl==0: break
                n = self.conn.recv(nl).hex()
                cl = int.from_bytes(self.conn.recv(4), 'big')
                c = self.conn.recv(cl).hex()
                print(f"\r[Alice]: {CryptoManager.decrypt(self.aes_key, n, c)}\nBob> ", end="")
        except Exception as e:
            print(f"\n[BOB RX] Connection closed: {e}")
            pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Bob Server - ECDH Receiver')
    parser.add_argument('testcase', type=str, help='Path to test case file (e.g., test_cases/35bit/case_1.txt)')
    args = parser.parse_args()
    
    BobServer(args.testcase).run()