#ifndef NUM_FORMAT_H_
#define NUM_FORMAT_H_

// str must be 26 + 3 + 1 + num decimals + 1 = 31+decimals bytes
// breakdown:
// strlen('18,446,744,073,709,551,615') = 26
// strlen(' GB') = 3
// strlen('.') = 1
// +1 for '\0'
char* bytes_to_str(unsigned long num, int decimals, char* str);

// Number to string using G to mean 10^9, M to mean 10^6 etc
char* num_to_str(double num, int decimals, char* str);


#endif /* NUM_FORMAT_H_ */
