from sys import argv
import sys


###USAGE###
#python3 calc_avg_depth file1.depth file2.depth



def getDepth(depthfile):
        '''
        Parses SGD features flat file
        Input: SGD_features.tab file
        Output: Avg depth
        '''

        lines=0
        totaldepth=0

        f = open(depthfile)
        for line in f:
                lines+=1
                row_data =line.split("\t")
                totaldepth+=int(row_data[2].strip())


        try:
                average_depth=totaldepth/lines
        except:
                average_depth=0
            
        return average_depth



if __name__ == '__main__':

        files = argv[1:]


        for each_file in files:

                depth=getDepth(each_file)
                print("depth: "+str(depth)+", 2xD: "+str(2*depth))
