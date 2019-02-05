#include <stdio.h>

void construct_input()
{
	//aux inputs
	//aux input W 64 words per instance
	//aux input abcdefgh 134 words per instance
	//total 198 words per instance with padding
	//n_leaf * 2 - 1 instances
	//total 198 * (n_leaf * 2 - 1) words in aux input
	//real input length = n_leaf * (512 / 32) words
	
	//last 256 bits final answer as aux
	
	//total inputs as bits 198 * 32 * (n_leaf * 2 - 1) + 512 * n_leaf = 13184 * n_leaf - 6080 bits
	
	
}

int main(int argc, char** argv)
{
	int n_leaf;
	sscanf("%d", &n_leaf);
	
	construct_input();
	
	return 0;
}
