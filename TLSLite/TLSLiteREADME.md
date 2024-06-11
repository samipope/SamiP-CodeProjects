# TLSLite Project
The TLSLite project is a custom implementation of a lightweight TLS (Transport Layer Security) protocol in Java, designed to secure communications between a client and a server. This project focuses on demonstrating the ability to establish secure connections, encrypt and decrypt messages, and ensure data integrity through digital signatures and MACs (Message Authentication Codes).

## Technologies used
- Java
- Cryptography Libraries (Java's built in libraries)
- Network Programming: Java sockets for TCP/IP communication

## Key Features
- Secure Socket Communication: Establishes secure socket connections between a client and a server using custom TLS-like handshakes.
- Encryption & Decryption: Implements AES encryption for securing messages transmitted between the client and the server.
- Certificate-Based Authentication: Uses X.509 certificates for server authentication to the client, ensuring that communications are with a trusted server.
- Diffie-Hellman Key Exchange: Implements the Diffie-Hellman protocol for secure key exchange, allowing both parties to derive a shared secret key.
- Digital Signatures & MACs: Ensures integrity and non-repudiation of messages using RSA digital signatures and HMACs.

 <img width="1142" alt="Screenshot 2024-06-10 at 9 51 44 PM" src="https://github.com/samipope/SamiP-CodeProjects/assets/142822253/9dad0107-bb4a-4792-bc55-4c9e3dd1cee7">


## How to run
1. Set up two terminals that can run concurrently
2. Enter the command "javac Server.java" in one terminal and "javac Client.java" in the other to compile both files
3. Start the server in one terminal using "java Server". Start the client in the other using "java Client".
   
