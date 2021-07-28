import os
import subprocess
import glob

def split_multiseg_file(filename):
    newFile = 0
    baseDir = os.path.dirname(filename)
    with open(filename) as f:
        for line in f.readlines():
            if line[0] == '>': # this is the start of a new feature, so make a new file                
                if newFile: # close old file if it's open
                    newFile.close()
                segName = line.strip()[5:-1] # this probably won't work in the general case
                print(segName)
                newFile = open(baseDir + '/' + segName + '.txt','w')
                continue
            else:
                newFile.write(line)
    return

def kml2gmt(filename):
    return

def resample_line(filename):
    # kind of a hacky way to use the subprocess module, but whatever
    cmd = 'gmt sample1d ' + filename + ' -T1m -AR > temp'
    subprocess.call(cmd,shell=True)
    subprocess.call('mv temp ' + filename,shell=True)
    return

if __name__ == '__main__':
    #regions = ['atlantic','indian','nazca']
    #for region in regions:
    #    split_multiseg_file('nazc-anta/nazca.txt')
    for file in glob.glob('nazc-anta/*R.txt'):
        resample_line(file)