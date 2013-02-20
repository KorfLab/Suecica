#!/usr/bin/env python3.2
import re, argparse
from sys import platform

def command_line():
	parser = argparse.ArgumentParser(description="DupsToClipboard.py opens a results file produced by DupHMM.py and copies all duplication \
									 locations to the clipboard in an IGV-friendly format. Coverage information is also included")
	parser.add_argument('-i', default="", type=str, help='DupHMM.py results filename.', metavar='ResultsInputFile')
	args = parser.parse_args()
	i_fn = args.i
	return(i_fn)
	
i_fn = command_line()
dups = []
pattern = r"^(\S+)\t(\d+)\t(\d+)\t\d+\t(\d+)\t(\S+)\t(\d+)\t(\d+)$"
recomp = re.compile(pattern)
with open(i_fn) as infile:
	for line in infile:
		match = recomp.match(line)
		if match:
			(chr, start, end, state, cov, TE_genes, TEs) = match.groups()
			if state == "2":
				outline = str(chr) + ":" + str(start) + "-" + str(end) + "\t(" + str(cov) + ")\t(" + str(TE_genes) + " TE genes)"
				if int(TE_genes) == 1: outline = outline[:-2] + ")"
				outline = outline + "\t(" + str(TEs) + " TEs)"
				if int(TE_genes) == 1: outline = outline[:-2] + ")"
				dups.append(outline)

if platform in ["win32","win64"]: # Windows
	cliptext = '\r\n'.join(dups)
	import win32clipboard
	win32clipboard.OpenClipboard()
	win32clipboard.EmptyClipboard()
	win32clipboard.SetClipboardData(win32clipboard.CF_UNICODETEXT,cliptext)
	win32clipboard.CloseClipboard()
else: # Linux and Mac
	cliptext = '\n'.join(dups)
	from subprocess import Popen, PIPE
	if "linux" in platform: # Linux
		xsel_proc = Popen(['xsel', '-pi'], stdin=PIPE)
		xsel_proc.communicate(bytearray(cliptext,"UTF-8"))
		xsel_proc = Popen(['xsel', '-bi'], stdin=PIPE)
		xsel_proc.communicate(bytearray(cliptext,"UTF-8"))
	elif platform in ["darwin", "os2", "os2emx"]: # Mac
		# Untested
		outtext = popen("pbcopy", "w")
		outtext.write(cliptext)
		outtext.close()
print(cliptext)