import os
os.system('./build.sh')
os.system('g++ -std=c++11 parser_sha_data_parallel.cpp -o psdp -O3')
for i in range(8):
	os.system('./psdp SHA256_64.pws SHA256_64_merkle_' + str(i + 1) + '_rdl.pws SHA256_64_merkle_' + str(i + 1) + '_circuit.txt SHA256_64_merkle_' + str(i + 1) + '_meta.txt')

os.system('make -C ../.. linear_gkr_zk')
os.system('cp ../../bin/main_zk .')
