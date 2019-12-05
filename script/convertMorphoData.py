# coding=UTF-8

# Converts a csv table with station, date, taxon and individuals / m2 to
# rows with station, date and all taxa

import sys

class station:
    """A station at a spec. date with all collected data"""
    def __init__(self, name, date):
        self.name = name
        self.dates = {date:{}}

    def addValue(self, date, taxon, value):
        if not date in self.dates:
            self.dates[date] = {}
        self.dates[date][taxon] = value


    def getDateData(self, date):
        return self.dates[date]
        

def main():


    # Read metadata
    
    stations = {}
    taxa = []
    
    md = open(sys.argv[1],"r")
    firstLine = True
    for line in md:
        if firstLine:
            firstLine = False
        else:
            items = line[:-1].split(",")
            st = items[0]
            date = items[1]
            taxon = items[2]
            value = items[3]

            if not st in stations:
                stations[st] = station(name = st, date=date)
            stations[st].addValue(date=date, taxon=taxon, value=value)

            if not taxon in taxa:
                taxa.append(taxon)

    md.close()
    
    
    # Print data

    print "\t".join(["Station", "Date"] + taxa)

    for sn in stations:
        st = stations[sn]
        for date in st.dates:
            pl = [st.name, date]
            dd = st.getDateData(date)

            for pn in taxa:
                if pn in dd:
                    pl.append(dd[pn])
                else:
                    pl.append("0")
            print "\t".join(str(f) for f in pl)
        
    
if __name__ == "__main__":

    main()
    
