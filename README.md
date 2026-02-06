Defelopment of a huffman encoder, main objective: being fast.
upto 64 bit on everything:
* huffmancodes: 2 to 64 bit, package algorithm ensures maximum huffman length, tested for 2 to 64 bit huffman lengths
* smybol set sizes: 1 to 63 bit symbol sizes. Huge symbolsize use a huge amount of memory, tested to 24 bit symbols
* symbol frequency range 2 to 64 bit
uses package algorithm when the standard huffman algorithm produces too long huffman codes
Some checking on the generated huffman codes is done.
Developed for 64 bit machines but can be easily configured for systems with less bits.
