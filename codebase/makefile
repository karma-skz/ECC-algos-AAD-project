PYTHON = python3
VENV = venv
ACTIVATE = source $(VENV)/bin/activate
TESTCASE = test_cases/35bit/case_5.txt

all: server

venv:
	$(PYTHON) -m venv $(VENV)
	$(ACTIVATE) && pip install -r requirements.txt

server: venv
	$(ACTIVATE) && $(PYTHON) web_interface/eve_server.py

alice: venv
	$(ACTIVATE) && $(PYTHON) web_interface/alice_server.py $(TESTCASE)

bob: venv
	$(ACTIVATE) && $(PYTHON) web_interface/bob_server.py $(TESTCASE)

clean:
	rm -rf web_interface/eve_temp.txt