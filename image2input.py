#!/usr/bin/env python

import sys, os, re, Image

def main(args):
	if len(args) < 3:
		print >>sys.stderr, "Usage: %s input_image output_text degree_of_connectivity=4 threashold=0" % os.path.basename(args[0])
		sys.exit(1)

	im = Image.open(args[1])
	f = file(args[2], "w")
	rows = im.size[0]
	cols = im.size[1]
	deg = 4
	th = 0
	if len(args) > 3: deg = int(args[3])
	if len(args) > 4: th = int(args[4])
	print >>f, im.size[0], deg, th
	for gray in im.getdata():
		if gray > 150:
			print >>f, 0,
		else:
			print>>f, 255,
	print >>f
	print ("done\n")
	f.close()
if __name__ == "__main__": main(sys.argv)
