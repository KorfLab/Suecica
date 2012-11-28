#! /usr/bin/env python3.2
import os, sys, math, re, threading, time; from shutil import copyfile

class Align(threading.Thread):
	def __init__(self, database, path, filename, pypath, locks):
		self.database = database
		self.path = path
		self.filename = filename
		self.pypath = pypath
		self.locks = locks
		threading.Thread.__init__(self)
	
	def run(self):
		global finishcount
		global totalthreads
		finished = 'F'
		file = ""
		while True:
			if finished == 'T': break
			counter = -1
			for lock in self.locks:
				if finished == 'T': break
				counter += 1
				if not lock.locked():
					with lock:
						file = str(self.path) + str(self.filename)
						self.path = "-p " + str(self.path)
						self.filename = "-f " + str(self.filename)
						print("\nStarting: " + str(file) + "\n")
						os.spawnl(os.P_WAIT,self.pypath, ('python', 'readcount_search.py', str(self.path), str(self.filename)))
						finished = 'T'
			time.sleep(5)
		finishcount += 1
		percentdone = round(finishcount / totalthreads * 100,0)
		print("\nFinished: " + str(file) + "\n" + str(percentdone) + "% done\n")

def queue(path):
	global finishcount
	global totalthreads
	filelist = []
	pypath = sys.executable
	database = "Data/t9_sorted.gff"
	dirlist = os.listdir(path)
	pattern = r'.*\.sam'
	recomp = re.compile(pattern)
	for filename in dirlist:
		match = recomp.search(filename)
		if match:
			filelist.append(filename)
	
	threadlist = []
	locks = []
	for i in range(0,5):
		if i < len(filelist):
			lock = threading.Lock()
			locks.append(lock)
	
	finishcount = 0
	totalthreads = len(filelist)
	
	for filename in filelist:
		p_thread = Align(database, path, filename, pypath, locks)
		p_thread.start()
		threadlist.append(p_thread)
	
	for p_thread in threadlist:
		p_thread.join()

start = time.time()

finishcount = 0
totalthreads = 0
path = "Data/SueAll/Sue-OurRESCAN/libsSue/ThalyrataAlignment/"
queue(path)
#path = "Data/SynNatLibraries/ProcessedF2s/2ndAlignment/Files430F2s/fqfiles/ThalyrataAlignment/"
#queue(path)
path = "Data/SynNatLibraries/ParentalData/Ler/lib1Ler/ThalyrataAlignment/"
queue(path)
path = "Data/SynNatLibraries/ParentalData/Ler/lib2Ler/ThalyrataAlignment/"
queue(path)
#path = "Data/SueAll/Sue-OurRESCAN/lib1Sue/ThalyrataAlignment/"
#queue(path)
#path = "Data/SueAll/Sue-OurRESCAN/lib2Sue/ThalyrataAlignment/"
#queue(path)
#path = "Data/CareAll/CareOurRESCAN/lib1Care1/ThalyrataAlignment/"
#queue(path)
#path = "Data/CareAll/CareOurRESCAN/lib2Care1/ThalyrataAlignment/"
#queue(path)

total = (time.time() - start) / 60
print("\nAlignments completed (" + str(total) + " minutes)")