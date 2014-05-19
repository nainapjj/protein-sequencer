import scipy
import csv

DATA_FILE = "parkinsons_updrs.data"

class Data:
    class DataColumn:
        def __init__(self, header):
            self.headerName = header
            self.dataRows = []
        
        def addData(self, dataPoint):
            self.dataRows.append(dataPoint)
    
    def __init__(self, dataFile):
        self.firstLine = ""
        self.otherLines = []
        self.numColumns = -1
        self.columnNames = []
        self.dataColumns = []
        self.__extractDataFile(dataFile)
    
    def setHeader(self, row):
        # Check to see if the number of columns changed.
        if (self.numColumns != -1 and self.numColumns != len(row)):
            print "Number of columns cannot be changed!"
        else:
            self.columnNames = row
            # If the number of columns is -1, then 
            # there shouldn't be data points yet.
            if (self.numColumns == -1):
                for column in row:
                    self.dataColumns.append(self.DataColumn(column)) 
            else:
                counter = 0
                for column in row:
                    self.dataColumns[counter].headerName = column
                    counter += 1
            
            self.numColumns = len(row)
    
    def addDataRow(self, row):
        # Check to see if the number of columns in the data row
        # is the same as the number of columns in the header
        if (self.numColumns != len(row)):
            print "Number of columns in data row does not match " + \
                "number columns in the header."
        else:
            counter = 0
            for dataPoint in row:
                self.dataColumns[counter].addData(dataPoint)
                counter += 1
    
    def printColumnByIndex(self, index):
        dColumn = self.getColumnByIndex(index)
        
        print self.dColumn.headerName
        print self.dColumn.dataRows
    
    def printColumnByName(self, name):
        dColumn = self.getColumnByName(name)
        
        if dColumn != None:
            print self.dColumn.headerName
            print self.dColumn.dataRows
    
    def getColumnByHeaderName(self, name):
        for column in self.dataColumns:
            if column.headerName == name:
                return column
        print "Column name not found."
        return None
    
    def getColumnByIndex(self, index):
        return self.dataColumns[index]
        
        
    def __extractDataFile(self, dataFile):
        with open(dataFile, 'r') as f:
            currentLineNumber = 0
            csvreader = csv.reader(f)
            
            for row in csvreader:
                currentLineNumber += 1
                currentLine = ", ".join(row)
            
                if (currentLineNumber == 1):
                    self.firstLine = currentLine
                    self.setHeader(row)
                    print "Found file with %d columns." % len(row)
                else:
                    self.otherLines.append(currentLine)
                    self.addDataRow(row)
            
            print "Done processing file. %d data rows found." % len(self.otherLines)
            
if __name__ == "__main__":
    p_data = Data()
    p_data.extractDataFile(DATA_FILE)
    #p_data.printColumnByIndex(2)