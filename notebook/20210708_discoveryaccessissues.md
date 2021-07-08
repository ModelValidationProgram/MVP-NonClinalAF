When Alan started working on his git branch in Discovery, he cloned that branch to Discovery. This prevented Katie from getting access.


First, Katie did this:
```
1. Go to the terminal and copy your public key file content.
cat ~/.ssh/id_rsa.pub
(If you don't have a key pair generated already, generate a new one using command "ssh-keygen" )
2. Go to the github website and login with your credentials.
3. In the upper-right corner of the github page, click your profile photo, then click Settings.
4. Click New SSH key or Add SSH key.
5. Paste your key into the "Key" field.
6. Click Add SSH key.
7. Test the connection by running the following command in the terminal.
ssh -T git@github.com 
```

Get message:
"Hi ModelValidationProgram/MVP-NonClinalAF! You've successfully authenticated, but GitHub does not provide shell access."

Then, Katie navigated to the repo on Discovery
```
cd /work/lotterhos/MVP-NonClinalAF
```
This seemed to work
