


import javax.crypto.BadPaddingException;
import javax.crypto.IllegalBlockSizeException;
import javax.crypto.NoSuchPaddingException;
import java.io.*;
import java.math.BigInteger;
import java.net.Socket;
import java.security.*;
import java.security.cert.Certificate;
import java.security.cert.CertificateException;
import java.util.Arrays;

public class Client {
    static BigInteger clientDHPrivateKey;
    static BigInteger sharedSecret;

public static void main(String[] Args) throws IOException, ClassNotFoundException, CertificateException, NoSuchAlgorithmException, SignatureException, InvalidKeyException, NoSuchProviderException, InvalidAlgorithmParameterException, NoSuchPaddingException, IllegalBlockSizeException, BadPaddingException {
    // -------------------- SETTING UP HANDSHAKE ----------------------

    //On each new connection from client to server, the server and client will pick random secret keys to use with Diffie-Hellman, which we'll call serverDHPriv for the server and clientDHPriv for the client. The public keys can be derieved: serverDHPub = g^serverDHPriv mod N and clientDHPub = g^clientDHPriv mod N.
    // client send a nonce to start
    BigInteger nonce =  new BigInteger(new SecureRandom().generateSeed(32)); //32 bytes for tls
    //make a ByteArrayOutputStream to keep the message history
    ByteArrayOutputStream messageHistory = new ByteArrayOutputStream();

    Socket clientSock =  new Socket("localhost",8080); //throws IOException
    ObjectOutputStream clientOutput = new ObjectOutputStream(clientSock.getOutputStream()); //make it an OBJECT output stream so you can send the nonce

    clientOutput.writeObject(nonce);
    messageHistory.write(nonce.toByteArray()); //write it to bytes
    System.out.println("Client sending the nonce!!");

    System.out.println("current working directory: " + System.getProperty("user.dir"));

    ObjectInputStream clientInput = new ObjectInputStream(clientSock.getInputStream());
    System.out.println("client reading in the first message from server....");
    Certificate serverCert = (Certificate) clientInput.readObject();
    BigInteger serverDHPublicKey = (BigInteger) clientInput.readObject();
    BigInteger signedServerDHPublicKey = (BigInteger) clientInput.readObject();

    Shared.verifyCertificate(serverCert);

    messageHistory.write(serverCert.toString().getBytes());
    messageHistory.write(serverDHPublicKey.toByteArray());
    messageHistory.write(signedServerDHPublicKey.toByteArray());

    clientDHPrivateKey = new BigInteger(new SecureRandom().generateSeed(32));
    sharedSecret = Shared.getSharedSecret(serverDHPublicKey,clientDHPrivateKey);
    System.out.println("client shared secret:" + sharedSecret);

    BigInteger clientDHPublicKey = Shared.getDHPublicKey(clientDHPrivateKey);


    File certFile = new File("../../CAcertificate.pem");
    System.out.println(certFile.getAbsolutePath());

    Certificate clientCert = Shared.getCertificate("../../CASignedClientCertificate.pem");
    PublicKey clientRSAPublicKey = clientCert.getPublicKey();
    BigInteger SIGNEDclientDHPublicKey = Shared.getSignedDHPublicKey("../../clientPrivateKey.der",clientDHPublicKey,clientRSAPublicKey);

    clientOutput.writeObject(clientCert);
    clientOutput.writeObject(clientDHPublicKey);
    clientOutput.writeObject(SIGNEDclientDHPublicKey);
    System.out.println("client sending certificate, public key and signed public key to server.....");

    messageHistory.write(clientCert.toString().getBytes());
    messageHistory.write(clientDHPublicKey.toByteArray());
    messageHistory.write(SIGNEDclientDHPublicKey.toByteArray());

    Shared.makeSecretKeys(nonce,sharedSecret);

    byte[] serverHandshakeFinishMsg = (byte[]) clientInput.readObject();
    Shared.checkForValidMAC(serverHandshakeFinishMsg, messageHistory.toByteArray(), Shared.serverMAC);

    messageHistory.write(serverHandshakeFinishMsg);
    System.out.println("client reading in serverHandshakeFinishMsg");

    byte[] clientHandshakeFinishMsg = Shared.macMessage(messageHistory.toByteArray(), Shared.clientMAC);
    clientOutput.writeObject(clientHandshakeFinishMsg);
    messageHistory.write(clientHandshakeFinishMsg);
    System.out.println("client sending handshake finish message");

    //----------------FINISHED HANDSHAKE-------------------------
    System.out.println("client handshake was successful!!");

    //--------------START MESSAGING ------------------------
    //read in the first message from the server side
    byte [] messagesFromServer = (byte[]) clientInput.readObject(); //read in the first obj
    String decryptedMessage = Shared.decrypt(messagesFromServer, Shared.serverEncryptionKey, Shared.serverInitVector, Shared.serverMAC);
    System.out.println("CLIENT-SIDE: Message from server (encrypted): " + Arrays.toString(messagesFromServer));
    System.out.println("CLIENT-SIDE: Message from server (decrypted): " + decryptedMessage);

    //then acknowledge that we go the message from the server
    String ACK = "ACK";
    byte[] encryptedACK = Shared.encrypt(ACK.getBytes(), Shared.clientEncryptionKey,Shared.clientInitVector, Shared.clientMAC);
    clientOutput.writeObject(encryptedACK); //send encrypted ACK
    System.out.println("CLIENT-SIDE: ACK message from me (decrypted): " + ACK);
    System.out.println("CLIENT-SIDE: ACK message from me (encrypted): " + Arrays.toString(encryptedACK));

    //then send a message that isn't the ack
    String ClientMessage = "Hey there server how you doin";
    byte[] encryptedClientMessage = Shared.encrypt(ClientMessage.getBytes(), Shared.clientEncryptionKey,Shared.clientInitVector,Shared.clientMAC);
    clientOutput.writeObject(encryptedClientMessage);
    System.out.println("CLIENT-SIDE: First real message from me (decrypted): " + ClientMessage);
    System.out.println("CLIENT-SIDE: First real message from me (encrypted): " + Arrays.toString(encryptedClientMessage));

    //receive the ack message that the server sent
    byte[] ACKFromServer =(byte[]) clientInput.readObject();
    String decryptedACKFromServer = Shared.decrypt(ACKFromServer,Shared.serverEncryptionKey,Shared.serverInitVector,Shared.serverMAC);
    System.out.println("CLIENT-SIDE: ACK message from Server (encrypted): " + Arrays.toString(ACKFromServer));
    System.out.println("CLIENT-SIDE: ACK message from Server (decrypted): " + decryptedACKFromServer);


if(!decryptedACKFromServer.equals("ACK")){
    throw new RuntimeException("CLIENT SIDE: ACK never received from the server");
}


} //end of main
} //end of class
