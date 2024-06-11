import javax.crypto.BadPaddingException;
import javax.crypto.IllegalBlockSizeException;
import javax.crypto.NoSuchPaddingException;
import java.math.BigInteger;
import java.net.ServerSocket;
import java.net.Socket;
import java.security.*;
import java.security.cert.Certificate;
import java.security.cert.CertificateException;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import java.util.Arrays;


public class Server {
    static BigInteger serverDHPrivateKey;
    static BigInteger sharedSecret;
    public static void main(String[] args) throws IOException, ClassNotFoundException, CertificateException, NoSuchAlgorithmException, SignatureException, InvalidKeyException, NoSuchProviderException, InvalidAlgorithmParameterException, NoSuchPaddingException, IllegalBlockSizeException, BadPaddingException {
        // -------------------- SETTING UP HANDSHAKE ----------------------
        ServerSocket serverSocket = new ServerSocket(8080);
        System.out.println("Server looking for connection from client....");
        Socket socket = serverSocket.accept();
        ObjectInputStream serverInput = new ObjectInputStream(socket.getInputStream()); //object to get nonce

        ByteArrayOutputStream messageHistory = new ByteArrayOutputStream();

        BigInteger nonce = (BigInteger) serverInput.readObject(); //throws classNotFoundException
        messageHistory.write(nonce.toByteArray());
        System.out.println("Server reading in nonce");
        //generate diffie-helman key
        serverDHPrivateKey = new BigInteger(new SecureRandom().generateSeed(32));
        //send the server certificate,
        ObjectOutputStream serverOutput = new ObjectOutputStream(socket.getOutputStream());

        BigInteger serverDHPublicKey = Shared.getDHPublicKey(serverDHPrivateKey);

        Certificate serverCert = Shared.getCertificate("../../CASignedServerCertificate.pem");

        PublicKey serverRSAPublicKey = serverCert.getPublicKey();
        BigInteger SIGNEDServerDHPublicKey = Shared.getSignedDHPublicKey("../../serverPrivateKey.der", serverDHPublicKey, serverRSAPublicKey);
        System.out.println("server sending second message including certifications and keys.....");

        serverOutput.writeObject(serverCert);
        serverOutput.writeObject(serverDHPublicKey);
        serverOutput.writeObject(SIGNEDServerDHPublicKey);

        messageHistory.write(serverCert.toString().getBytes());
        messageHistory.write(serverDHPublicKey.toByteArray());
        messageHistory.write(SIGNEDServerDHPublicKey.toByteArray());

        System.out.println("server is reading in client's certificate, public key, and signed public key");
        Certificate clientCert = (Certificate) serverInput.readObject();
        BigInteger clientDHPublicKey = (BigInteger) serverInput.readObject();
        BigInteger SIGNEDClientDHPublicKey = (BigInteger) serverInput.readObject();

        Shared.verifyCertificate(clientCert);

        messageHistory.write(clientCert.toString().getBytes());
        messageHistory.write(clientDHPublicKey.toByteArray());
        messageHistory.write(SIGNEDClientDHPublicKey.toByteArray());

        sharedSecret = Shared.getSharedSecret(clientDHPublicKey,serverDHPrivateKey);
        System.out.println("server shared secret:" + sharedSecret);

        Shared.makeSecretKeys(nonce,sharedSecret);

        byte[] serverHandshakeFinishMsg = Shared.macMessage(messageHistory.toByteArray(),Shared.serverMAC);
        serverOutput.writeObject(serverHandshakeFinishMsg);
        messageHistory.write(serverHandshakeFinishMsg);
        System.out.println("server sending handshakefinishMSG.....");

        byte [] clientHandshakeFinishMsg = (byte[]) serverInput.readObject();
        System.out.println("server reading client handshakefinishMSG");

        Shared.checkForValidMAC(clientHandshakeFinishMsg, messageHistory.toByteArray(), Shared.clientMAC);
        messageHistory.write(clientHandshakeFinishMsg);
        // -------------------- FINISHED HANDSHAKE ---------------------------------
        System.out.println("server handshake was successful!!!");

        //--------------- START MESSAGING --------------------------

        //server sends the first message
        String firstMessage = "Hello client, let us start talking and stuff";
        byte[] encryptedFirstMessage = Shared.encrypt(firstMessage.getBytes(),Shared.serverEncryptionKey,Shared.serverInitVector,Shared.serverMAC);
        serverOutput.writeObject(encryptedFirstMessage);
        System.out.println("SERVER-SIDE: First message from me (decrypted): " + firstMessage);
        System.out.println("SERVER-SIDE: First message from me (encrypted): " + Arrays.toString(encryptedFirstMessage));

        //receive the ACK from client
        byte[] ACKFromClient = (byte[]) serverInput.readObject();
        String decryptedACKFromClient = Shared.decrypt(ACKFromClient, Shared.clientEncryptionKey, Shared.clientInitVector,Shared.clientMAC);
        System.out.println("SERVER-SIDE: ACK from client (encrypted): " + Arrays.toString(ACKFromClient));
        System.out.println("SERVER-SIDE: ACK from client (decrypted): " + decryptedACKFromClient);

        if(!decryptedACKFromClient.equals("ACK")){
            throw new RuntimeException("SERVER SIDE: ACK never received from the client");
        } //check to make sure the ACK is working properly

        //read in the message from the client (should == "Hey there server how you doin")
        byte[] messageFromClient = (byte[]) serverInput.readObject();
        String decryptedMessageFromClient = Shared.decrypt(messageFromClient, Shared.clientEncryptionKey, Shared.clientInitVector,Shared.clientMAC);
        System.out.println("SERVER-SIDE: Message from client (encrypted): " + Arrays.toString(messageFromClient));
        System.out.println("SERVER-SIDE: Message from client (decrypted): " + decryptedMessageFromClient);

        //send ACK back
        String ACK = "ACK";
        byte[] encryptedACK = Shared.encrypt(ACK.getBytes(),Shared.serverEncryptionKey,Shared.serverInitVector,Shared.serverMAC);
        serverOutput.writeObject(encryptedACK); //send it over
        System.out.println("SERVER-SIDE: ACK message from me (decrypted): " + ACK);
        System.out.println("SERVER-SIDE: ACK message from me (encrypted): " + Arrays.toString(encryptedACK));



    }
}
