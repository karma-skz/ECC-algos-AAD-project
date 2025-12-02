## How to setup :

```bash

# Create virtual environment
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt


# Start the server
python eve_server.py


# Webpage
Visit http://127.0.0.1:5500
Click "Start Capture"

# Open another terminal and start Alice and Bob clients
python3 web_interface/alice_server.py <testcase_path>
python3 web_interface/bob_server.py <testcase_path>

#Select the handshake packet and Click on "FullAttack" to crack the ECDH key

# Chat, then select the new packet and Decrypt messages!

```