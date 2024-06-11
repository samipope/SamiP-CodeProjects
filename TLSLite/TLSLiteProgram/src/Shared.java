
import javax.crypto.*;
import javax.crypto.spec.IvParameterSpec;
import javax.crypto.spec.SecretKeySpec;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.security.*;
import java.security.cert.Certificate;
import java.security.cert.CertificateException;
import java.security.cert.CertificateFactory;
import java.security.spec.InvalidKeySpecException;
import java.security.spec.PKCS8EncodedKeySpec;
import java.util.Arrays;


public class Shared {

    //I used a safeprime number that is from the Kivinen and Kojo standards. It is 1536-bit and comes from this website: https://datatracker.ietf.org/doc/html/rfc3526#section-2
    static private final BigInteger n = new BigInteger("FFFFFFFFFFFFFFFFC90FDAA22168C234C4C6628B80DC1CD129024E088A67CC74020BBEA63B139B22514A08798E3404DD" +
            "EF9519B3CD3A431B302B0A6DF25F14374FE1356D6D51C245" +
            "E485B576625E7EC6F44C42E9A637ED6B0BFF5CB6F406B7ED" +
            "EE386BFB5A899FA5AE9F24117C4B1FE649286651ECE45B3D" +
            "C2007CB8A163BF0598DA48361C55D39A69163FA8FD24CF5F" +
            "83655D23DCA3AD961C62F356208552BB9ED529077096966D" +
            "670C354E4ABC9804F1746C08CA237327FFFFFFFFFFFFFFFF", 16);
    static public BigInteger g = new BigInteger("2"); // "The generator is: 2" is from the website
    static Certificate CACert = Shared.getCertificate("../../CAcertificate.pem");

    //SESSION KEYS
    static byte[] serverEncryptionKey; //session specific key for server
    static byte[] clientEncryptionKey; //session specific key for client
    static byte[] serverMAC; //server message authentication codes
    static byte[] clientMAC; // client message authentication codes
    static byte[] serverInitVector; //used by server when encrypting and sending data to client. ensures it is unique
    static byte[] clientInitVector; //used by the client when encrypting and sending data to server. ensures it is unique


    /**
     * generates a certificate
     * @param filepath --> passing in the path of the client signed CA file
     * @return the certificate made by the certificate factory built in functionality
     */
    static public Certificate getCertificate(String filepath) {
        try {FileInputStream fileInput = new FileInputStream(filepath);
            CertificateFactory certFact = CertificateFactory.getInstance("X.509"); //"x.509" is a parameter that tells the CertFact that it is TLS
            return certFact.generateCertificate(fileInput);
        }catch (FileNotFoundException | CertificateException e){
            throw new RuntimeException(e);
        }
    }

    /**
     * DIFFIE-HELLMAN --> calculates the public key from private key
     * @param DHprivateKey private key
     * @return public key
     */
    public static BigInteger getDHPublicKey(BigInteger DHprivateKey) {
        return g.modPow(DHprivateKey, n);
    }

    /**
     * DIFFIE-HELLMAN --> signs a diffie-helman public key using RSA and verifies signiature
     * @param privateKeyDer_file
     * @param DHpublicKey
     * @param RSAPublicKey
     * @return
     */
    static public BigInteger getSignedDHPublicKey(String privateKeyDer_file, BigInteger DHpublicKey, PublicKey RSAPublicKey){
        try (FileInputStream privateKeyDer_file_input = new FileInputStream(privateKeyDer_file)) {
            byte[] privateKeyBytes = privateKeyDer_file_input.readAllBytes();
            PKCS8EncodedKeySpec encodedPrivateKey = new PKCS8EncodedKeySpec(privateKeyBytes);
            KeyFactory keyFactory = KeyFactory.getInstance("RSA");
            PrivateKey privateRSAKey = keyFactory.generatePrivate(encodedPrivateKey);
            Signature signature = Signature.getInstance("SHA256withRSA"); //creates signature object using 256-bit hash value and RSA
            signature.initSign(privateRSAKey);
            signature.update(DHpublicKey.toByteArray());
            byte[] signedPublicKey = signature.sign();
            signature.initVerify(RSAPublicKey);
            signature.update(DHpublicKey.toByteArray());
            return new BigInteger(signedPublicKey);
        } catch (IOException | NoSuchAlgorithmException | InvalidKeySpecException | InvalidKeyException | SignatureException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * verifies the certificate against the CA's public key
     * @param certificate
     * @throws CertificateException
     * @throws NoSuchAlgorithmException
     * @throws SignatureException
     * @throws InvalidKeyException
     * @throws NoSuchProviderException
     */
    public static void verifyCertificate(Certificate certificate) throws CertificateException, NoSuchAlgorithmException, SignatureException, InvalidKeyException, NoSuchProviderException {
        PublicKey _CApublicKey = Shared.CACert.getPublicKey();
        certificate.verify(_CApublicKey);
    }


    /**
     * calculates the SHARED secret using the DIFFIE-HELLMAN key exchange
     * @param publicDHKey
     * @param privateDHKey
     * @return
     */
    public static BigInteger getSharedSecret(BigInteger publicDHKey, BigInteger privateDHKey) {
        return publicDHKey.modPow(privateDHKey, n);
    }

    /**
     * expands a master key into a derived key --> uses HMAC-based key derived function
     * @param masterKey
     * @param tag
     * @return byte array of derived key that is a fixed length
     */
    private static byte[] hdkfExpand(byte[] masterKey, String tag) {
        try {
            Mac HMAC = Mac.getInstance("HmacSHA256");
            SecretKeySpec secretKeySpec = new SecretKeySpec(masterKey, "HmacSHA256");
            HMAC.init(secretKeySpec);
            HMAC.update(tag.getBytes(StandardCharsets.UTF_8));
            HMAC.update((byte) 0x01);
            return Arrays.copyOf(HMAC.doFinal(), 16);
        } catch (NoSuchAlgorithmException | InvalidKeyException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * makes a session key for encryption and mac
     * @param clientNonce
     * @param sharedSecretFromDiffieHellman
     */
    public static void makeSecretKeys(BigInteger clientNonce, BigInteger sharedSecretFromDiffieHellman) {
        byte[] prk = hdkfExpand(clientNonce.toByteArray(), sharedSecretFromDiffieHellman.toString());
        serverEncryptionKey = hdkfExpand(prk, "server encrypt");
        clientEncryptionKey = hdkfExpand(serverEncryptionKey, "client encrypt");
        serverMAC = hdkfExpand(clientEncryptionKey, "server MAC");
        clientMAC = hdkfExpand(serverMAC, "client MAC");
        serverInitVector = hdkfExpand(clientMAC, "server IV");
        clientInitVector = hdkfExpand(serverInitVector, "client IV");
    }

    /**
     * makes a MAC for a given message using the key u pass in
     * @param message
     * @param macKey
     * @return byte array
     */
    public static byte[] macMessage(byte[] message, byte[] macKey) throws NoSuchAlgorithmException, InvalidKeyException, IOException {
        Mac HMAC = Mac.getInstance("HmacSHA256");
        SecretKeySpec secretKeySpec = new SecretKeySpec(macKey, "HmacSHA256");
        HMAC.init(secretKeySpec);
        HMAC.update(message);
        return HMAC.doFinal();
    }

    /**
     *
     * @param sendersMacMsg
     * @param myMessageHistory
     * @param macKey
     * @throws NoSuchAlgorithmException
     * @throws IOException
     * @throws InvalidKeyException
     */
    public static void checkForValidMAC(byte[] sendersMacMsg, byte[] myMessageHistory, byte[] macKey) throws NoSuchAlgorithmException, IOException, InvalidKeyException {
        byte[] myMacMsg = Shared.macMessage(myMessageHistory, macKey);
        if (!Arrays.equals(myMacMsg, sendersMacMsg)){
            throw new RuntimeException("Message MACS mismatch");
        }

    }

    /**
     * creates and initializes link cipher instance for encryption or decryption operations
     * @param isEncryptCipher
     * @param encryptKey
     * @param initializationVector
     * @return new cipher
     *
     */
    public static Cipher createCipher( Boolean isEncryptCipher, byte[] encryptKey, byte[] initializationVector) throws NoSuchPaddingException, NoSuchAlgorithmException, InvalidAlgorithmParameterException, InvalidKeyException {
        Cipher cipher = Cipher.getInstance("AES/CBC/PKCS5Padding");
        SecretKeySpec secretKeySpec = new SecretKeySpec(encryptKey, "AES");
        IvParameterSpec ivParameterSpec = new IvParameterSpec(initializationVector);
        cipher.init(isEncryptCipher ? Cipher.ENCRYPT_MODE : Cipher.DECRYPT_MODE, secretKeySpec, ivParameterSpec);
        return cipher;
    }


    /**
     * encrypts a message for either the server or client depending on which encryption key, init vector and mac key are sent
     * @param message
     * @param encryptKey
     * @param initializationVector
     * @param macKey
     * @return byte array that has encrypted your string
     */
    public static byte[] encrypt(byte[] message, byte[] encryptKey, byte[] initializationVector, byte[] macKey) throws IOException, NoSuchAlgorithmException, InvalidKeyException, NoSuchPaddingException, InvalidAlgorithmParameterException, IllegalBlockSizeException, BadPaddingException {
        ByteArrayOutputStream encryptedMessage = new ByteArrayOutputStream();
        encryptedMessage.write(message);
        byte[] macMsg = macMessage(message, macKey);
        encryptedMessage.write(macMsg);
        Cipher cipher = createCipher(true, encryptKey, initializationVector);
        return cipher.doFinal(encryptedMessage.toByteArray());

    }

    /**
     * decrypts a message for either the server or client depending on which encryption key, init vector and mac key are sent
     * @param cipherText
     * @param encryptKey
     * @param initializationVector
     * @param macKey
     * @return String that has been decrypted
     */
    public static String decrypt(byte[] cipherText, byte[] encryptKey, byte[] initializationVector, byte[] macKey) throws InvalidAlgorithmParameterException, NoSuchPaddingException, NoSuchAlgorithmException, InvalidKeyException, IllegalBlockSizeException, BadPaddingException {
        Cipher cipher = createCipher(false, encryptKey, initializationVector);
        byte[] plainText = cipher.doFinal(cipherText);
        byte[] decryptedMsg = Arrays.copyOf(plainText, plainText.length - 32);
        return new String(decryptedMsg, StandardCharsets.UTF_8);
    }
}







