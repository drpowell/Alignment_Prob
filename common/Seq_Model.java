
package common;

import java.io.*;
/**
   The interface for any left to right model of sequence 
*/
public interface Seq_Model extends Serializable {
    /** Calculate what it would cost to encode a character
	@see #update
	@param a The character to encode
	@param i The character <b>a</b> is the <b>i</b><sup>th</sup> character of the sequence
	@return The length to encode the character <b>a</b>
    */
    double encodeLen(char a, int i);

    /** Update the internal model for encoding a character.  Like {@link #encodeLen} but actually updates
	the internal state of the model as required.
	@see #encodeLen
	@param a The character to encode
	@param i The character <b>a</b> is the <b>i</b><sup>th</sup> character of the sequence
	@return The length to encode the character <b>a</b>
    */
    double update(char a, int i);
}

