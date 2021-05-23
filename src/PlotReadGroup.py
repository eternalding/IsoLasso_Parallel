import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import argparse
import time
import os
import matplotlib._color_data as mcd
import tqdm


class RG():
    def __init__(self):
        self.Name           = ""
        self.ChrName        = ""
        self.Range          = []
        self.ExonBoundary   = []
        self.ExonCoverage   = []
        self.SGTypes        = []
        self.TypeCount      = []
        self.ValidReads     = 0
        self.Orientation    = "." 
        self.TypeDir        = []
        self.CvgStats       = []


def ReadFile(fp):

    line = fp.readline()[:-1]
    while line:
        if(line[0] == "*"):
            line = fp.readline()[:-1]
            continue
        else:
            break

    ReadGroups = []
    NewRG = RG()

    while(line):
        if line.find("ReadGroup")!=-1:
            if(NewRG.Name!=""):
                ReadGroups.append(NewRG)
            NewRG = RG()
            _,NewRG.Name = line.split(" ")
            line = fp.readline()[:-1]
        elif line.find("Boundary")!=-1:
            _,NewRG.ChrName,start,end,NewRG.Orientation = line.split(" ")
            NewRG.Range.append(int(start))
            NewRG.Range.append(int(end))
            line = fp.readline()[:-1]
        elif line.find("Total valid reads:")!=-1:
            NewRG.ValidReads = int(line[line.find(":")+1:])
            line = fp.readline()[:-1]
            Numexon = int(line[line.find(":")+1:])
            line = fp.readline()[:-1]

            linecount = 0
            while(linecount<Numexon):
                Start,End,ExonLen,Coverage,Cvgstats = line.split(" ")
                NewRG.ExonBoundary.append((int(Start),int(End)))
                NewRG.ExonCoverage.append(int(Coverage))
                NewRG.CvgStats.append(Cvgstats)
                line = fp.readline()[:-1]
                linecount+=1

        elif line.find("SGTypes")!=-1:
            NumSGTypes = int(line[line.find(":")+1:])
            Typeindex = 0
            line = fp.readline()[:-1]
            while(Typeindex<NumSGTypes):
                TypeInfo = line.split(" ")
                NewRG.SGTypes.append(TypeInfo[:-2])
                NewRG.TypeCount.append(int(TypeInfo[-2]))
                NewRG.TypeDir.append(TypeInfo[-1])
                line = fp.readline()[:-1]
                Typeindex+=1
        else:
            line = fp.readline()[:-1]    


    ReadGroups.append(NewRG)
    return ReadGroups

def parabola(point1, point2,max):
    x = np.linspace(point1[0], point2[0], 100)

    interval = 0.5*(point2[0]-point1[0])

    inc = np.linspace(-interval , interval , 100)
    
    y = (-max/interval**2)*inc**2+max
    return (x,y+point1[1])

parser =argparse.ArgumentParser()
parser.add_argument("-i","--input" ,help="The input file of RGStats.txt")
parser.add_argument("-n","--numreads" ,help="The minimum reads count for a readgroup to be displayed, default:1000",default=1000,type=int)
parser.add_argument("-r","--RGperPic" ,help="The number of readgroups displayed per picture.",default=4,type=int)
parser.add_argument("-m","--minTypeCnt" ,help="The minimum required typecount for a type to be displayed.",default=100,type=int)
args = parser.parse_args()

input_filename = args.input[args.input.rfind("/")+1:]

local_time = time.localtime()
if not os.path.isdir("../result/"):
    os.mkdir("../result/")

result_path = "../result/"+str(local_time.tm_year)+"."+str(local_time.tm_mon)+"."+str(local_time.tm_mday)+" "+str(local_time.tm_hour)+":"+str(local_time.tm_min)+":"+str(local_time.tm_sec)+" "+str(input_filename[:-4])+"/"

os.mkdir(result_path)

fp = open(args.input,'r')

RGs = ReadFile(fp)

index = 0

fig = plt.figure(figsize=(120, 60),dpi=80)
plt.tight_layout()

Subplot_index = 1 
Pic_index = 0

font_dict = {'family':'monospace',
            'color':'black',
            'size':20}

cmap =plt.get_cmap("tab20")
fig.suptitle('IsoLasso ReadGroup Analysis', fontsize=80)
plt.subplots_adjust(hspace=1)


for RG in tqdm.tqdm(RGs):
    
    if RG.ValidReads<int(args.numreads):
        continue

    axis = fig.add_subplot(args.RGperPic,1,Subplot_index)
    axis.set_title("Readgroup "+RG.Name,fontsize=40) 
    axis.set_ylabel(RG.ChrName,fontsize=40)
    
    for tick in axis.xaxis.get_major_ticks():
        tick.label.set_fontsize(36)

    axis.ticklabel_format(style='plain')

    axis.set_xlim((int(RG.ExonBoundary[0][0]),int(RG.ExonBoundary[-1][1])))
    axis.axhline(0.5,color='black',ls='-')

    exon_color = mcd.XKCD_COLORS["xkcd:maroon"].upper()

    CvgSum = sum(RG.ExonCoverage)
    exon_pos = []


    for exon_index in range(len(RG.ExonBoundary)):

        alpha = int(RG.ExonCoverage[exon_index]*10/CvgSum)/10 
        if alpha==0:
            alpha=0.1

        exon_patch =axis.add_patch(
                                patches.Rectangle(
                                (int(RG.ExonBoundary[exon_index][0]),0.45),
                                int(RG.ExonBoundary[exon_index][1])-int(RG.ExonBoundary[exon_index][0]),
                                0.1,
                                edgecolor = 'black',
                                facecolor = exon_color,
                                fill=True,
                                zorder =100,
                                label='exon',
                                alpha= 1
                            ))

        axis.text(RG.ExonBoundary[exon_index][0],0.6,RG.ExonBoundary[exon_index][0],horizontalalignment='center',verticalalignment='center',fontdict=font_dict)
        axis.text(RG.ExonBoundary[exon_index][1],0.4,RG.ExonBoundary[exon_index][1],horizontalalignment='center',verticalalignment='center',fontdict=font_dict)


        exon_pos.append(RG.ExonBoundary[exon_index][0])
        exon_pos.append(RG.ExonBoundary[exon_index][1])
        
    TypeCount_Sum = sum(RG.TypeCount)
    direction = 1

    for idx, Type in enumerate(RG.SGTypes):
    
        if RG.TypeCount[idx]>args.minTypeCnt:
            Exons = [i for i, x in enumerate(Type) if x==str(1)]

            alpha = int(RG.TypeCount[idx]*100/TypeCount_Sum)*2 
            if alpha==0:
                alpha=1

            if len(Exons)==1:
                X,y = parabola((0.25*int(RG.ExonBoundary[Exons[0]][0])+0.75*int(RG.ExonBoundary[Exons[0]][1]),0.5),(0.75*int(RG.ExonBoundary[Exons[0]][0])+0.25*int(RG.ExonBoundary[Exons[0]][1]),0.5),direction*alpha/100)
                axis.plot(X,y,color=cmap(idx%20)[:3],linewidth=alpha)
                axis.text(X[49],y[40],RG.TypeCount[idx],horizontalalignment='center',fontdict=font_dict)

            else:
                for exon_idx in range(len(Exons)):
                    if exon_idx==len(Exons)-1:
                        break
                    else:
                        X,y = parabola((RG.ExonBoundary[Exons[exon_idx]][1],0.5),(RG.ExonBoundary[Exons[exon_idx+1]][0],0.5),direction*alpha/100)
                        axis.plot(X,y,color=cmap(idx%20)[:3],linewidth=alpha)
                        direction*=-1

                    axis.text(X[49],y[40],RG.TypeCount[idx],horizontalalignment='center',fontdict=font_dict)
            

    axis.legend(handles=[exon_patch],loc="upper right",prop={'size': 30})

    Subplot_index+=1
    if Subplot_index>args.RGperPic:
        Subplot_index = 1
        plt.savefig(result_path+"ReadGroup_"+str(Pic_index)+".jpg")
        Pic_index+=1
        plt.clf()

#fig.colorbar(im, ax=axis, shrink=0.6)

plt.savefig(result_path+"ReadGroup_"+str(Pic_index)+".jpg")


