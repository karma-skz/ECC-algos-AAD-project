import os
import hashlib
import json
from cryptography.hazmat.primitives.ciphers.aead import AESGCM

class CryptoManager:
    """
    Handles Real-World AES-256-GCM Encryption.
    """
    
    @staticmethod
    def derive_key(shared_point_x: int) -> bytes:
        """Derive a 32-byte AES key from the Shared Secret X-coord."""
        # SHA-256 KDF (Standard practice)
        return hashlib.sha256(str(shared_point_x).encode()).digest()

    @staticmethod
    def encrypt(key: bytes, plaintext: str) -> dict:
        """Encrypts text -> AES-GCM Ciphertext."""
        try:
            aesgcm = AESGCM(key)
            nonce = os.urandom(12)
            data = plaintext.encode()
            ct = aesgcm.encrypt(nonce, data, None)
            return {
                'nonce': nonce.hex(),
                'ciphertext': ct.hex()
            }
        except Exception as e:
            print(f"Crypto Error: {e}")
            return None #type:ignore

    @staticmethod
    def decrypt(key: bytes, nonce_hex: str, ct_hex: str) -> str:
        """Decrypts Hex -> Plaintext."""
        try:
            aesgcm = AESGCM(key)
            nonce = bytes.fromhex(nonce_hex)
            ct = bytes.fromhex(ct_hex)
            return aesgcm.decrypt(nonce, ct, None).decode()
        except:
            return "[DECRYPTION FAILED]"