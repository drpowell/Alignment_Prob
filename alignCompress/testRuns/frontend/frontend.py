#!/usr/bin/python

import os,sys,re,random
import bisect,string, types
import signal,cPickle
import time
from qt import *

# Global vars
data_fd=0
data_str=''
numArch=0
numPop=0
startTime=0
upto=0
remainingTime=0
type=0
mainForm=None

def makeUI(filename):
    s = os.popen("pyuic %s" % filename)
    code = s.read()
    s.close()
    return code

def openFile(fname):
    fd = os.open(fname, os.O_RDONLY)
    return fd

def readFile(fd):
    #print "In readFile"
    str = ''
    while 1:
        s = os.read(fd, 10000000)
        if len(s)==0: break
        str += s
    #print "Leaving readFile"
    return str

def pollFile():
    global data_fd
    str = readFile(data_fd)
    if len(str):
        #print ':%s:' % str[:-1]
        update(str)

def readLogFile(fname):
    #print "read "+fname
    global numArch, numPop, startTime
    try:
        f = open(fname)
    except IOError:
        print "Unable to find log file '"+fname+"'.  Guessing stuff..."
        numArch=10
        numPop=200
        startTime=time.time()
        return
    s = f.read()
    strs = re.split(r"];", s, 2)

    numArch = 0
    numPop = 0
    for which in [0,1]:
        for l in re.split(r"}", strs[which]):
            l = re.sub(r"^(.|\n)*?{", "{", l)
            if l[0] != '{': continue
            l = re.sub(r"=>", r":", l)
            l = 'info = ' + l + '}'
            exec(l)
            sub_pos = info['S_LEN']
            e_pos = info['S_LEN'] + info['SUB_LEN']
            seq = info['SEQ'][0:sub_pos-1] + ' ### ' + info['SEQ'][sub_pos:e_pos-1] + ' ### ' + info['SEQ'][e_pos:]

            if which==0:
                mainForm.Sequences.append("Parent %03d: %s"%(numArch,seq))
                numArch = numArch + 1
            else:
                mainForm.Sequences.append("Child %03d: %s"%(numPop,seq))
                numPop = numPop + 1
        
    m = re.search(r"STDERR:(\d+):", s)
    startTime = int(m.group(1))

def updateProgress(str):
    global upto, numPop, numArch, remainingTime
    m = re.compile(r".*PRSS: s1=(\d+) s2=(\d+)", re.DOTALL).match(str)
    if not m: return
    a = int(m.group(1))
    p = int(m.group(2))
    if (a*numPop+p+1 > upto):
        upto = a*numPop+p+1
        mainForm.ProgressBar1.setProgress(upto)

        eta = (time.time()-startTime)*(numPop*numArch)/(upto) + startTime
        remainingTime = int(eta - time.time())

        txt = "Done %d of %d\n" % (upto, numArch*numPop)
        txt += "Started at %s\n"%(time.strftime("%d/%m/%Y %H:%M:%S",time.localtime(startTime)))
        txt += "Now is %s\n"%( time.strftime("%d/%m/%Y %H:%M:%S",time.localtime(time.time())))
        txt += "ETA %s\n"%( time.strftime("%d/%m/%Y %H:%M:%S",time.localtime(eta)))
        mainForm.ProgStatus.setText(txt)

        displayRemainingTime()

def displayRemainingTime():
    mainForm.LCD.display("%02d:%02d:%02d" %
                         (remainingTime/3600, (remainingTime/60)%60, remainingTime%60))

def countDown():
    global remainingTime
    if (remainingTime>0):
        remainingTime -= 1
    displayRemainingTime()

def update(str):
    global data_str
    data_str += str
    updateProgress(str)
    left = mainForm.Status.xOffset()
    top = mainForm.Status.topCell()
    lastVis = (mainForm.Status.lastRowVisible() >= mainForm.Status.numRows()-1)
    mainForm.Status.append(str[:-1])
    if lastVis:
        mainForm.Status.setTopCell(mainForm.Status.numRows())
    else:
        mainForm.Status.setTopCell(top)
    mainForm.Status.setXOffset(left)
    


def plotProg():
    if type==0: return
    pid = os.fork()
    if pid:
        #parent
        return
    prog = (os.path.dirname(resFile) or '.') + '/populationPlot' + ('.pl','2.pl')[type-1]
#    os.close(1)                     # No stdout for sub-process
    os.execv(prog, [prog, resFile])
    sys.exit(1)._exit()


def main(argv):
    global mainForm, data_fd, type, resFile
    
    if len(argv) != 2:
        print "Usage: %s <res.pop file>" % argv[0]
        sys.exit(1)

    resFile = argv[1]
    #resFile = '../res.population2.20020917.txt'

    rundir = os.path.dirname(argv[0])
    if not rundir:
        rundir='.'
    exec(makeUI(rundir + '/' + 'gui.ui'))

    app = QApplication(sys.argv)
    mainForm = MainForm()

    data_fd  = openFile(resFile)
    str = readFile(data_fd)

    m = re.compile(r"^#pid=(\d+)", re.MULTILINE).search(str)

    if not m:
        print "Unable to find 'pid=' line"
        sys.exit(1)

    pid = m.group(1)

    if re.search(r"Blend Model", str):
        type = 2
    else:
        type = 1

    d = os.path.dirname(resFile) or "."
    readLogFile( d + '/' + 'popLog' + ( ('.','2.')[type-1]) + pid )

    

    mainForm.ProgressBar1.setTotalSteps(numArch * numPop)

    update(str)
    mainForm.Status.setLeftCell(0)

    mainForm.setCaption(resFile + '   Type=%d' % type)
    mainForm.setIconText('Output Viewer')
    QObject.connect(mainForm.PlotButton,  SIGNAL( 'clicked()' ), plotProg)

    mp3_poll = QTimer()
    QObject.connect(mp3_poll, SIGNAL( 'timeout()' ), pollFile)
    mp3_poll.start(500)

    countDownTimer = QTimer()
    QObject.connect(countDownTimer, SIGNAL( 'timeout()' ), countDown)
    countDownTimer.start(1000)

    #mainForm.move(0,0)
    app.setMainWidget(mainForm)
    mainForm.show()
    app.exec_loop()


main(sys.argv)







