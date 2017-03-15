#ifndef ALTRUISM_H
#define ALTRUISM_H

#include <iostream>
#include <cmath>
using namespace std;
// main altruism functions
// always run in the following order only

int main(){

	
	for (int i=0; i<=10; ++i){
		float x = 0+i*50/10;	
//		float p_disperse = 1/(1+exp(x - 15));
		float p_disperse = 1/(1+exp(10*(x - 15)));
		cout << x << "\t" << p_disperse << endl;
	}
	
	return 0;
}


#endif
