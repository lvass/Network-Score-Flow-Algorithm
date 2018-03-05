import copy
import re

#function to read in network data
def readingCSV (interactionsFile, scoresFile):
    print('reading file start')
    inAdj = []
    outAdj = []
    sign = []
    fileObj = open(interactionsFile, 'r')
    linelist = fileObj.readlines()
    del linelist[0]
    pattern = ','
    for interaction in linelist: #reading in csv file
        col = re.split(pattern, interaction)
        newCol = []
        for item in col:
            newItem = item.replace('\r\n', '')
            newCol.append(newItem)
        inAdj.append(newCol[0])
        outAdj.append(newCol[1])
        sign.append(int(newCol[2]))
    
    scoreDict= {} #creating dictionaries to store the raw scores
    blankDict = {} #dictionary with all nodes equal to score 0 to store results in
    scoresData = open(scoresFile, 'r')
    linelist = scoresData.readlines()
    del linelist[0]
    pattern = ','
    for score in linelist:
        col = re.split(pattern, score) #splitting into columns 
        newCol = []
        for item in col:
            newItem = item.replace('\r\n', '')
            newCol.append(newItem)
        scoreDict[str(newCol[0])] = float(newCol[1]) #values equal to raw scores from csv file
        blankDict[str(newCol[0])] = int(0) #all values equal to 0
    scoreDict.update({'P53':10}) #use to change p53 raw score value
    print('reading file end')
    return(inAdj, outAdj, sign, scoreDict, blankDict)

#function to levelize the network
def levelCalc(inAdj, outAdj, cycleValue):
    print('Algorithm 1 start')
    l0 = [] #looking for 'start' nodes, those 
    for item in inAdj:
        try: 
            outAdj.index(item)
        except ValueError:
            try: 
                l0.index(item)
            except ValueError:
                l0.append(item)
    start = l0
    l0 = [l0]
    levels = []
    visitedList= []
    counter = 0
    repeat = True
    while repeat == True:
        i = 0
        thisLevel = []
        for startNode in start:
            cycle = False
            if len(start) == 0:
                cycle = True
            position = [index for index, value in enumerate(inAdj) if value == startNode] #looking for the node in inAdj
            dupPosition = [index for index, value in enumerate(visitedList) if value == startNode] #looks for duplicates to identify cycles/multiple in degrees
            if len(dupPosition) > cycleValue: #limits the number of cycles and max in degree of a node
                cycle = True
            if cycle == False:
                for place in position:
                    thisLevel.append(outAdj[place]) #adding node to the level
                    visitedList.append(outAdj[place]) #counts visits
            if cycle == True:
                pass
            i+=1
        if len(thisLevel) ==0:
            repeat = False    #stops loop when no more levels exist
        else:
            levels.append(thisLevel) #adding level to the list of levels
            start = thisLevel
            counter +=1
    
    l0.extend(levels) #adding the start levels onto levels list
    print('Algorithm 1 end')
    return(l0)

#function to calculate the scores in a step-wise manner. Stops when end point score thresholds are reached
def scoreCalc(inAdj, outAdj, levels, scoreDict, sign, blankDict, endNodes, thresholds):
    print('Algorithm 2 start')
    resultsDict = copy.deepcopy(blankDict) #all values 0
    visitedDict = copy.deepcopy(blankDict) #all values 0
    for node in levels[0]:
        resultsDict[node]+= scoreDict[node] #initiating start level scores
    stepCounter = 0
    interactionCount = 0 #prevents endless cycles, limits to cycle limit value
    i = 0
    thresholdReached = 0
    while (i < len(levels)):
        for node in levels[i]:
            instances = [index for index, value in enumerate(inAdj) if value == node] #looking for node in inAdj list
            totalAdjScores = 0
            interactionFromList = []
            interactionToList =[]
            signList =[]
            for x in instances:
                interaction = outAdj[x]
                totalAdjScores+= scoreDict.get(interaction) #calculating total raw score of all out adjacent nodes of the node 'x'
                interactionFromList.append(inAdj[x]) #creating a nodes from which the score is travelling from
                interactionToList.append(outAdj[x]) #creating a nodes from which the score is travelling to
                signList = signList + [sign[x]]
            z = 0
            while z < len(interactionFromList): 
                stepCounter+=1
                mySign = signList[z]
                raw = float(scoreDict.get(interactionToList[z]))
                output = float(resultsDict.get(interactionFromList[z]))
                totalAdj = float(totalAdjScores)
                current = float(resultsDict.get(interactionToList[z]))
                try:
                    myScore = float(mySign * raw * (output/totalAdj)) #calculating score as described by Isik 2012
                except ZeroDivisionError:
                    z+=1
                    break
                if visitedDict.get(interactionToList[z]) > 0:
                    (resultsDict[interactionToList[z]])+= myScore #if not first visit, don't re-add raw score
                else:
                    interactionCount+=1
                    (resultsDict[interactionToList[z]])+= myScore + raw #if first visit, also add raw score
                    visitedDict[interactionToList[z]] +=1
                thresholdOfNode = (thresholds.get(interactionToList[z])) #checking if in end node list
                if thresholdOfNode != None:
                    if thresholdOfNode <= resultsDict[interactionToList[z]]-1: #checking if threshold has been met
                        thresholdReached+=1
                        if thresholdReached == 1: #first end point threshold met adds 1 to thresholdReached
                            firstEnd = str(interactionToList[z]) #saving result
                            firstStep = stepCounter
                            thresholds.pop(firstEnd) #deleting from end point list
                            break
                if thresholdReached == 2: # second end point threshold met makes thresholdReached equal 2
                    secondEnd = str(interactionToList[z]) #saving results
                    secondStep = stepCounter
                    break #breaking out of loops and ending function
                z+=1
            if thresholdReached == 2:
                break #breaking out of loops and ending function
        i+=1
        if thresholdReached == 2:
            break #breaking out of loops and ending function
    return (resultsDict, firstEnd, firstStep, secondEnd, secondStep)



#Main program

data = (readingCSV('network_data.csv', 'dataset1.csv'))

my_endNodes = ['CELLULAR SENESCENCE', 'APOPTOSIS'] #identify end points
my_thresholds = {'CELLULAR SENESCENCE':(100000), 'APOPTOSIS':(100000)} #set parameters

my_levels = levelCalc(data[0], data[1], 156) #set cycle limit

finalScores = scoreCalc(data[0], data[1], my_levels, data[3], data[2], data[4], my_endNodes, my_thresholds)


print '{0} threshold score reached after {1} interactions'.format(finalScores[1], finalScores[2])
print '{0} threshold score reached after {1} interactions'.format(finalScores[3], finalScores[4])
print('end')



