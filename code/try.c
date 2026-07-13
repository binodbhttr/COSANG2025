#include <stdio.h>

int main(void) {

	unsigned int pid = 234567;
	long long mostbound;

	mostbound = pid;
	
	printf("%lld %u\n", mostbound, pid);
	
	if (mostbound == pid)
		printf("OK\n");

	return 0;

}
