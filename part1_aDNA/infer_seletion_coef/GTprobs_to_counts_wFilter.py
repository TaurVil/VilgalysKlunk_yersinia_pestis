'''
This script parses all .genolik files for given set of loci {LOCI; exons, gwas, neutral} and populations {POP; london, denmark}. Samples in pre, during, post are considered to occur 16, 8, and 1 generations ago for UK, and 18, 9, 1 for denmark. 

Within each panel, individuals with >50% missing GT calls will be removed & reported.

This script returns expectations of counts (sum of count*GTprobs ), which are floats.

usage: python GTprobs_to_counts_wFilter.py POP LOCI

'''
import sys, re, time
import numpy as np
int_regex = re.compile(r'[0-9]+') # regex expression to grep all integers


def genolik2expct(line):
	'''
	Convert triplets of genolik to expected count / mean count
	Return a list of counts and a list of total indexed by individual
	'''
	# sanity check
	try:
		assert len(line) % 3 == 0
	except Exception as e:
		print(f'len(line) = {len(line)}. {e}.')
		sys.exit()
	# get number of ind.
	sampNum = len(line) / 3
	# going through (AA, Aa, aa) triplets
	count = []; total = []
	for i in range(int(sampNum)):
		# skip if no data
		if '-1.' in line[i*3:(i*3 + 3)]:
			count.append(0)
			total.append(0)
			continue
		lik_AA, lik_Aa, lik_aa = map(float, line[i*3:(i*3 + 3)] )
		count.append(2.*lik_AA + lik_Aa)
		total.append(2)

		# sanity check
		try:
			assert np.isclose(sum((lik_AA, lik_Aa, lik_aa)) , 1., atol = 1e-5)
		except Exception as e:
			print(f'sum(lik_AA, lik_Aa, lik_aa) = {sum((lik_AA, lik_Aa, lik_aa))}.')
			#sys.exit()
			continue

	if sum(total) == 0 and set(line) != {'-1.'}:
		print(count, total, line, lik_AA, lik_Aa, lik_aa)
		#sys.exit()
	return(count, total, int(sampNum))


def readGenolik_asMatrice(filename):
	'''Output a dictionary of GT matrix, {chr: [X matrix, N matrix, pos_list]}. x is the major allele for each locus'''
	Freq_by_chr = {}
	with open(filename , 'r') as geno:
		for l in geno:
			l = l.strip().split(' ')
			# get chromosome index
			chrom = int(int_regex.findall(l[0])[0])
			# get position
			pos = int(l[1])
			x, n, sampNum = genolik2expct(l[2:])
			if chrom -1 not in Freq_by_chr:
				Freq_by_chr[chrom -1] = [ np.empty((0,sampNum)), np.empty((0,sampNum)), [] ]
			Freq_by_chr[chrom - 1][0] = np.vstack( (Freq_by_chr[chrom - 1][0], np.array(x)) ) 
			Freq_by_chr[chrom - 1][1] = np.vstack( (Freq_by_chr[chrom - 1][1], np.array(n)) ) 
			Freq_by_chr[chrom - 1][2].append(pos)
	# make sure file is closed
	geno.close()
	return(Freq_by_chr)


def filter_missing(Freq_by_chr):
	'''Report and remove individuals with >50% GT calls genome-wide'''
	sampNum = Freq_by_chr[0][1].shape[1]
	#Filtered_Freqs = []
	# tally each chromosome
	total_loci = 0
	missing_loci = np.empty((0, sampNum))
	for c in range(22):
		#chrom = c+1
		N_Matrix = Freq_by_chr[c][1]
		# sanity check
		assert N_Matrix.shape[1] == sampNum
		total_loci += N_Matrix.shape[0] #nrow
		missing_loci = np.vstack((missing_loci, np.sum( (N_Matrix == 0), axis = 0)) ) #sum up all rows by each col
	print(f'Total number of locus: {total_loci}')
	# sum up all chromosome
	missing = np.sum( missing_loci, axis = 0 )
	# get the ids
	filtered_num = np.sum(missing >= 0.5*total_loci)
	filtered_ids = np.array(np.where(missing >= 0.5*total_loci)).reshape((filtered_num))
	# generate filtered Freqs
	ids_to_keep = np.array(np.where(missing < 0.5*total_loci)).reshape((sampNum - filtered_num))
	#print(np.shape(ids_to_keep))
	Filtered_Freqs = {}
	for c in range(22):
		X, N, pos_list = Freq_by_chr[c]
		numSites = X.shape[0]
		x_counts = np.sum(X[:,ids_to_keep], axis = 1).reshape((numSites))
		n_counts = np.sum(N[:,ids_to_keep], axis = 1).reshape((numSites))

		Filtered_Freqs[c] = [ x_counts, n_counts, pos_list ]
		'''try:
		except Exception as e:
			print(np.sum(X[:,ids_to_keep], axis = 1))
			print(X[:,ids_to_keep])
			print(X[:,ids_to_keep].shape)
			print(e)
			sys.exit()'''
	return(Filtered_Freqs, filtered_ids)


def write_all(Trajs, filtered_ids, outname):
	'''Write all counts to file & record the IDs of removed individuals/samples'''
	# simplify stuff
	def _formatIndices(ids):
		return( ', '.join( map(str, ids) ) )

	outfile = open(outname, 'w')
	# write sampling time and header
	if 'london' in outname:
		outfile.write(f'## Removed GT calls from {len(filtered_ids[0])} t1 individual: {_formatIndices(filtered_ids[0])}\n## Removed GT calls from {len(filtered_ids[1])} t2 individual {_formatIndices(filtered_ids[1])}\n## Removed GT calls from {len(filtered_ids[2])} t3 individual {_formatIndices(filtered_ids[2])}\n## Removed sites either with no samples across time or with MAF == 0.\n## Recording all sites otherwise\n##SampTimes.gen.ago: 16, 8, 1\nchr\tposition\tx1\tn1\tx2\tn2\tx3\tn3\n')

	elif 'denmark' in outname:
		outfile.write(f'## Removed GT calls from {len(filtered_ids[0])} t1 individual: {_formatIndices(filtered_ids[0])}\n## Removed GT calls from {len(filtered_ids[1])} t2 individual {_formatIndices(filtered_ids[1])}\n## Removed GT calls from {len(filtered_ids[2])} t3 individual {_formatIndices(filtered_ids[2])}\n## Removed sites either with no samples across time or with MAF == 0.\n## Recording all sites otherwise\n##SampTimes.gen.ago: 18, 9, 1\nchr\tposition\tx1\tn1\tx2\tn2\tx3\tn3\n')
	# start writing
	skipped = {}
	for c in range(1, 23):
		skipped[c] = []

		# retrieve the list of positions
		pos_list = Trajs[0][c-1][2]
		# sanity check
		try:
			assert set(pos_list) == set(Trajs[1][c-1][2])
			assert set(pos_list) == set(Trajs[2][c-1][2])
			assert set(Trajs[2][c-1][2]) == set(Trajs[1][c-1][2])
		except:
			print('Loci recorded at three time points do not match.')
			sys.exit()
		# no need to sort positions
		for i in range(len(pos_list)):
			pos = pos_list[i]
			# skip sites with no samples
			N = [ np.sum(Trajs[t][c-1][1][i]) for t in range(3) ]
			try:
				total_N = np.sum(N)
			except Exception as e:
				print(f'chrom {c}, pos {pos}\nN = {N}')
				print(N)
				print(e)
				sys.exit()
			if total_N == 0:
				skipped[c].append(pos)
				continue
			# skip the sites with no counts
			X = [ np.sum(Trajs[t][c-1][0][i]) for t in range(3)]
			total_X = np.sum(X)
			if total_X == total_N or total_X == 0:
				skipped[c].append(pos)
				#continue
			else:
				oline = '%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\n' % (c, pos, X[0], N[0], X[1], N[1], X[2], N[2])
				outfile.write(oline)		
		print(f'Removed {len(skipped[c])} non-segregating sites from chromosome {c}: {skipped[c]}')
		# finished
	outfile.close()



def main():
	import os
	path = os.getcwd()
	POP, LOCI = sys.argv[1:]

	print(f'{time.ctime()}. Start reading inputs for {POP} population on {LOCI} set.')
	Allele_trajs = [] ; Filtered_IDs = []
	for TIME in ('pre','during','post'):
		infile = os.path.join(path, "genoliks", f"genolik.{LOCI}_{POP}_{TIME}.genolik")
		Freqs = readGenolik_asMatrice(infile)
		Filtered_Freqs, filtered_ids = filter_missing(Freqs)
		Allele_trajs.append(Freqs)
		Filtered_IDs.append(filtered_ids)
		print(filtered_ids)

	# write output
	#outpre = f'{path}sample_counts/{POP}Samples_{LOCI}_filtered_WG.txt'
	outfile = os.path.join(path, f"{POP}Samples_{LOCI}_filtered_WG.txt")
	print(f'{time.ctime()}. Start writing output to {outfile}.')
	write_all(Allele_trajs, Filtered_IDs, outfile)
	print(f'{time.ctime()} Done.\n\n')

if __name__ == '__main__':
	main()