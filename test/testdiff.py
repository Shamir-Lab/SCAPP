

f = open('test_out/test.confident_cycs.fasta')
line = f.readlines()[0]
assert 'RNODE_1_length_4088_cov_433.' in line, "Test Failed!"
f.close()
print("Test passed")
