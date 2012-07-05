#! /usr/bin/env python3.2
import os, sys, math, re, threading, time; from shutil import copyfile

class Align(threading.Thread):
	def __init__(self, database, file, pypath, locks):
		self.database = database
		self.file = file
		self.pypath = pypath
		self.locks = locks
		threading.Thread.__init__(self)
	
	def run(self):
		#usage: bwa-samtools.py [-q trim quality] [-pb path to bwa] [-ps path to samtools] database_file.fa reads.fq
		global finishcount
		global totalthreads
		finished = 'F'
		while True:
			if finished == 'T': break
			counter = -1
			for lock in self.locks:
				if finished == 'T': break
				counter += 1
				if not lock.locked():
					with lock:
						newdatabase = str(self.database[:-3]) + "-" + str(counter) + ".fa"
						print("\nStarting:\n" + str(self.file) + "\n(" + str(newdatabase) + ")\n")
						os.spawnv(os.P_WAIT,self.pypath, ('python', 'bwa-samtools.py', str(newdatabase), str(self.file)))
						os.remove(str(newdatabase) + ".amb")
						os.remove(str(newdatabase) + ".ann")
						os.remove(str(newdatabase) + ".bwt")
						os.remove(str(newdatabase) + ".fai")
						os.remove(str(newdatabase) + ".pac")
						os.remove(str(newdatabase) + ".rbwt")
						os.remove(str(newdatabase) + ".rpac")
						os.remove(str(newdatabase) + ".rsa")
						os.remove(str(newdatabase) + ".sa")
						finished = 'T'
			time.sleep(5)
		finishcount += 1
		percentdone = round(finishcount / totalthreads * 100,1)
		print("\nFinished: " + str(self.file) + "\n" + str(percentdone) + "% done\n")

start = time.time()

filelist = []
pypath = sys.executable
database = "Data/Thalyrata-w-breaks.fa"
#database = "../at9_genome.fasta"
#path = "Data/SynNatLibraries/ProcessedF2s/2ndAlignment/Files430F2s/fqfiles/"
#path = "Data/SueAll/Sue-OurRESCAN/lib2Sue/"
#path = "Data/SueAll/Sue-OurRESCAN/lib1Sue/"
#path = "Data/SynNatLibraries/ParentalData/Ler/lib2Ler/"
#path = "Data/SynNatLibraries/ParentalData/Ler/lib1Ler/"
#path = "Data/CareAll/CareOurRESCAN/lib2Care1/"
#path = "Data/CareAll/CareOurRESCAN/lib1Care1/"
#Following list has 3 Sue, 3 Ler FASTQ files
#filelist = ["Data/SueAll/Sue-OurRESCAN/lib2Sue/libSue1.fq","Data/SueAll/Sue-OurRESCAN/lib1Sue/libSue1F.fq","Data/SueAll/Sue-OurRESCAN/lib1Sue/libSue1R.fq","Data/SynNatLibraries/ParentalData/Ler/lib1Ler/libLer1F.fq","Data/SynNatLibraries/ParentalData/Ler/lib1Ler/libLer1R.fq", "Data/SynNatLibraries/ParentalData/Ler/lib2Ler/libLer2.fq"]
filelist = ["Data/CareAll/CareOurRESCAN/lib2Care1/lib2Care1.fq","Data/CareAll/CareOurRESCAN/lib1Care1/libCare1F.fq","Data/CareAll/CareOurRESCAN/lib1Care1/libCare1R.fq"]
num_threads = 10
#dirlist = os.listdir(path)#
pattern = r'.*\.fq'
recomp = re.compile(pattern)
#for filename in dirlist:#
#	match = recomp.search(filename)#
#	if match:#
#		filelist.append(filename)#

threadlist = []
locks = []
for i in range(0,num_threads):
	if i < len(filelist):
		lock = threading.Lock()
		locks.append(lock)
		newdatabase = database[:-3] + "-" + str(i) + ".fa"
		if not os.path.exists(newdatabase):
			copyfile(database,newdatabase)

finishcount = 0
totalthreads = len(filelist)

for filename in filelist:
	#file = str(path) + filename
	file = filename
	p_thread = Align(database, file, pypath, locks)
	p_thread.start()
	threadlist.append(p_thread)

for p_thread in threadlist:
	p_thread.join()

for i in range(0,num_threads):
	if i < len(filelist):
		newdatabase = database[:-3] + "-" + str(i) + ".fa"
		os.remove(str(newdatabase))

total = (time.time() - start) / 60
print("\nAlignments completed (" + str(total) + " minutes)")