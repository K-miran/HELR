#ifndef SCHEME_SCHEME_H_
#define SCHEME_SCHEME_H_

#include "SecKey.h"
#include "Cipher.h"
#include "Message.h"
#include "PubKey.h"
#include "SchemeAux.h"


using namespace std;
using namespace NTL;
class Scheme {
private:
public:
    Params& params;
    SecKey& secretKey;
    PubKey& publicKey;
    SchemeAux& aux;
    
    //-----------------------------------------
    
    Scheme(Params& params, SecKey& secretKey, PubKey& publicKey, SchemeAux& schemeaux) : params(params), secretKey(secretKey), publicKey(publicKey), aux(schemeaux) {};
    //-----------------------------------------
    
    ZZ getqi(long& modulus);
    ZZ getPqi(long& modulus);
    long getLogqi(long& modulus);
    
    //-----------------------------------------
    
    void trueValue(CZZ& m, ZZ& qi);
    void trueValue(ZZ& m, ZZ& qi);
    
    //-----------------------------------------
    
    void rlweInstance(ZZX& bx, ZZX& ax, ZZ& qi);
    void rlweInstance(ZZX& bx, ZZX& ax);
    
    //-----------------------------------------
    
    CZZ* groupidx(CZZ*& vals, long slots);
    CZZ* groupidx(CZZ& val);
    CZZ* degroupidx(CZZ*& vals, long dslots);
    
    //-----------------------------------------
    
    Message encode(CZZ*& vals, long size);
    CZZ* decode(Message& msg);
    
    //-----------------------------------------
    
    Cipher encrypt(Message& msg, long modulus = 0);
    
    Message decrypt(Cipher& cipher);
    
    //-----------------------------------------
    
    Cipher encryptSlots(CZZ*& vals, long slots, long modulus = 0);
    CZZ* decryptSlots(Cipher& cipher);
    
    
    Cipher encrypt_small(Message& msg, long modulus = 0);  // perform modulus-switching before adding to the message
    Cipher encryptSlots_small(CZZ*& vals, long slots, long modulus = 0);
    //-----------------------------------------
    
    Cipher encryptSlots(CZZ& val, long modulus = 0);
    CZZ decryptFullSingle(Cipher& cipher);
    
    //-----------------------------------------
    
    Cipher* encryptFullSingleArray(CZZ*& vals, long size);
    CZZ* decryptFullSingleArray(Cipher*& ciphers, long size);
    
    //-----------------------------------------
    
    Cipher add(Cipher& cipher1, Cipher& cipher2);
    void addAndEqual(Cipher& cipher1, Cipher& cipher2);
    
    //-----------------------------------------
    
    Cipher addConst(Cipher& cipher, ZZ& cnst);
    void addConstAndEqual(Cipher& cipher, ZZ& cnst);
    
    //-----------------------------------------
    
    Cipher sub(Cipher& cipher1, Cipher& cipher2);
    void subAndEqual(Cipher& cipher1, Cipher& cipher2);
    
    Cipher subConst(Cipher& cipher, ZZ& cnst);
    void subConstAndEqual(Cipher& cipher, ZZ& cnst);  // cipher - cnst
    void subTransConstAndEqual(ZZ& cnst, Cipher& cipher);  // cnst - cipher
    
    //-----------------------------------------
    
    Cipher mult(Cipher& cipher1, Cipher& cipher2);
    void multAndEqual(Cipher& cipher1, Cipher& cipher2);
    
    //-----------------------------------------
    
    Cipher square(Cipher& cipher);
    void squareAndEqual(Cipher& cipher);
    
    //-----------------------------------------
    
    Cipher multByConst(Cipher& cipher, ZZ& cnst);
    void multByConstAndEqual(Cipher& cipher, ZZ& cnst);
    
    //-----------------------------------------
    
    Cipher multByMonomial(Cipher& cipher, const long& degree);
    void multByMonomialAndEqual(Cipher& cipher, const long& degree);
    
    //-----------------------------------------
    
    Cipher leftShift(Cipher& cipher, long& bits);
    void leftShiftAndEqual(Cipher& cipher, long& bits);
    void doubleAndEqual(Cipher& cipher);
    
    //-----------------------------------------
    
    //Cipher modSwitch(Cipher& cipher, long newLevel);
    //Cipher modSwitchOne(Cipher& cipher);
    //void modSwitchAndEqual(Cipher& cipher, long newLevel);
    //void modSwitchOneAndEqual(Cipher& cipher);
    
    
    Cipher modSwitch(Cipher& cipher, long scalebits);
    void modSwitchAndEqual(Cipher& cipher, long scalebits);
    
    //-----------------------------------------
    
    //Cipher modEmbed(Cipher& cipher, long newLevel);
    //Cipher modEmbed(Cipher& cipher);
    //void modEmbedAndEqual(Cipher& cipher, long newLevel);
    //void modEmbedOneAndEqual(Cipher& cipher);
    
    
    Cipher modEmbed(Cipher& cipher, long& newModulus);
    void modEmbedAndEqual(Cipher& cipher, long& newModulus);
    //-----------------------------------------
    
    Cipher multAndModSwitchOne(Cipher& cipher1, Cipher& cipher2);
    void multModSwitchOneAndEqual(Cipher& cipher1, Cipher& cipher2);
    
    //-----------------------------------------
    
    Cipher rotate2(Cipher& cipher, long& logSteps);
    void rotate2AndEqual(Cipher& cipher, long& logSteps);
    Cipher rotate(Cipher& cipher, long& steps);
    void rotateAndEqual(Cipher& cipher, long& steps);
    
    //-----------------------------------------
    
};


#endif /* SCHEME_SCHEME_H_ */
