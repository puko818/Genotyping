__author__ = 'Ivan Vogel'
# naming structure:
# for gDNA: "sampleNR_gDNA.Gtype"
# for trios: "sampleNR_trioNR [PB1|PB2|egg].Gtype"
#
#Name;Chr;Position


import sys
import csv
import re
import collections
import numpy as np
from tabulate import tabulate
from collections import Counter
import itertools



codes={
    'AA': 3,
    'AB': 2,
    'BB': 0,
    'NC': -1,
    'NONHET': -9,
    'AA/AB':4,
    'AA/AA':3,
    'AA/BB':6,
    'BB/BB':0,
    'BB/AB':4,
    'BB/AA':9
}


names={
    4: "HET",#AA/AB,BB/AB - red
    2: "HET",#code for gDNA - red
   -9: "NONHET",#if gDNA==NONHET - white
    3: "MATCH", #AA/AA, - orange
    0: "MATCH", #BB/BB - orange
    6: "MISM", #AA/BB - green
    9: "MISM", #"BB/AA" - green
    -1: "NOCALL",#gray
}


class Columnnames:
    snpname ='Name'
    chromosome='Chr'
    position='Position'
    frompos=3



'''
class Samples(object):
    def __init__(self):
        self._datasets=[]
    def add_rec(self,rec):
        self._datasets.append(rec)
'''


class ChromosomeCollection(object):
    ref = None
    gdna=None
    trios=None

    def __init__(self,gdna,reference,trios):
        self._col=dict()
        ChromosomeCollection.trios=trios
        ChromosomeCollection.gdna=gdna
        ChromosomeCollection.ref=reference

    class Chromosome:
        def __init__(self,chr_id):
            self._samples=dict()
            self.ID=chr_id
            self._hetpos =None


        def hetpos(self):
            if (not self._hetpos): self._hetpos=[1 if i==codes["AB"] else 0 for i in self._samples[ChromosomeCollection.gdna]]
            return self._hetpos

        def statistics(self):
            pass


        '''
        def __getitem__(self,item):
            return self._samples[item]
        '''

        def __getitem__(self, item):
            #return [call for call in self._samples[item] for hetpos in self._hetpos if hetpos==1]
           return [call for (call,hetpos) in zip(self._samples[item],self._hetpos) if hetpos==1]

        def __setitem__(self,sampleid,call):
            if sampleid not in self._samples.keys(): self._samples[sampleid]=[]
            self._samples[sampleid].append(call)

        def dumptolist(self): return [pos for sample in self._samples.values() for pos in sample]

        def dumptofile(self):
            #print tabulate([self._samples[k] for k in self._samples.keys()], headers=self._samples.keys())
            with open('testfile.txt', 'a') as the_file:
                the_file.write(tabulate([self._samples[k] for k in self._samples.keys()], headers=[str(i) for i in  range(0,len(self._samples.keys()))],tablefmt="plain"))

            #print header
            #print("\t".join(self._samples.keys()))
            #for k in self._samples.keys():
            #    print("\t".join([]))

        def _notref(self,s): return s!=ChromosomeCollection.ref

        def _notgdna(self,s): return s!=ChromosomeCollection.gdna

        def resolve_parents(self,call,gdna,ref):
            retval=codes['NONHET']
            if gdna==1:
                if ref==codes['AA']:
                  if call==codes['AB']:
                      retval=codes['AA/AB']

                  elif call==codes['AA']:
                      retval=codes['AA/AA']

                  elif call==codes['BB']:
                      retval=codes['AA/BB']

                elif ref==codes['BB']:
                  if call==codes['AB']:
                      retval=codes['BB/AB']
                  elif call==codes['AA']:
                      retval=codes['BB/AA']
                  elif call==codes['BB']:
                      retval=codes['BB/BB']
                if call==codes['NC']:
                    retval=codes['NC']
            return retval



        def phase(self):
            output=[]
            #for k in self._samples:
            for k in [el for trio in ChromosomeCollection.trios for el in trio]:
                if self._notgdna(k) or self._notref(k):
                    output.append((k,[self.resolve_parents(call,gdna,ref) for (call,gdna,ref) in zip(self._samples[k],self.hetpos(),self._samples[ChromosomeCollection.ref])]))
                    #output.append((k,[call for (call,gdna,ref) in zip(self._samples[k],self.hetpos(),self._samples[ChromosomeCollection.ref])]))
            return output

            #for body in [self._samples[k] for k in self._samples if self._notref(k) and self._notgdna(k)]



        def create_minimal_gdna_bed(self):
           startpos = 0
           #[1 if i==codes["AB"] else 0 for i in self._samples[ChromosomeCollection.gdna]]
           for k, g  in itertools.groupby(self._samples[ChromosomeCollection.gdna]):
           ###!!!!!call len(list(g)) only once! every call probably next() method
             currentrange=len(list(g))
             code=codes["NONHET"]
             if k==codes["AB"]:
               code=k
             print "%s\t%d\t%d\t%d\t%s"%(ChromosomeCollection.gdna,int(startpos),int(startpos+currentrange),int(code),names[code])
             startpos+=currentrange



        '''
        def getNCs(self,sampleid):
            count=0
            #return [i if i==-1 else i=0  for i in self._samples[sampleid]]
            return self._samples[sampleid].count(-1)
        '''



    def __str__(self): return ','.join(sorted([chr for chr in self._col.keys()]))

    '''
    def __setitem__(self, key, value):
        if chromosomeID not in self._col.keys():
           self._col[chromosomeID]=self.Chromosome(chromosomeID)
    '''
    def __getitem__(self, chromosomeID):
       if chromosomeID not in self._col.keys():
           self._col[chromosomeID]=self.Chromosome(chromosomeID)
       return self._col[chromosomeID]
    '''
    def __setitem__(self, chromosomeID, value):

        else:
            with self._col[chromosomeID] as chr:
                chr[]
    '''



class Call:
    Code={'AA': 1, 'BB':2, 'AB': 3, 'NC':-1}

    def __init__(self):
        self._call='NC'

    def get_call_code(self):
        pass

    def resolve_call(self):
        pass


Rec = collections.namedtuple('Sample', 'ColumnName Sample TrioID Type ColPos')



class MetaData(object):
    def __init__(self,headers):
      self._dataset=[]
      #Rec = collections.namedtuple('ColumnName', 'Sample TrioID Type ColPos')
      for idx, val in enumerate(headers):
        m1 = re.match(r"(?P<sampleid>\d+)_(?P<trioid>\d+)\s(?P<type>\bPB1\b|\bPB2\b|\begg\b).GType", val)
        #check if it's trio
        m2=None
        if m1:
          s,t,ty=m1.group('sampleid','trioid','type')
          self._dataset.append(Rec(ColumnName=val,Sample=s,TrioID=t,Type=ty,ColPos=idx))
        else:
          m2=re.match(r"(?P<sampleid>\d+)_(?P<type>gDNA|OVO|PB1).GType", val)
          if m2:
            s2,ty2=m2.group('sampleid','type')
            self._dataset.append(Rec(ColumnName=val,Sample=s2,TrioID=-1,Type=ty2,ColPos=idx))


    def get_sample(self,sampleid):
        return [d for d in self._dataset if int(d.Sample)== int(sampleid)]

    def get_gdna_by_sample(self,sampleid):
        return [':'.join((d.ColumnName,str(d.ColPos))) for d in self._dataset if int(d.Sample)== int(sampleid) and d.Type=="gDNA"][0]

    def get_trios_by_sample(self,sampleid):
        """

        :rtype : list of trio tuples of concat. ColumnName+ColumnPos
        """
        l= [':'.join((d.ColumnName,str(d.ColPos))) for d in self._dataset if int(d.Sample)== int(sampleid) and d.Type in ("PB1","PB2","egg")]
        return [s for s in zip(*[iter(l)]*3)]



    def get_reference(self,sampleid):
        return [':'.join((d.ColumnName,str(d.ColPos))) for d in self._dataset if int(d.Sample)== int(sampleid) and d.Type in ("PB2","egg")][0]


    def dump(self):
        print "%s\t%s\t%s\t%s\t%s" % ("Dataset","Sample","Trio","Type","Column Position in Input File")
        for d in self._dataset:
            print "%s\t%s\t%s\t%s\t%s" %  (d.ColumnName,d.Sample,d.TrioID,d.Type,d.ColPos)



def correct_keys(reader):
    """Correct reader for unique keys
    :param reader: CSV reader
    :param columns: columns with potential redundancy
    """
    for rara in reader:
        pairs={}
        for i,k in enumerate(reader.fieldnames):
            if i<Columnnames.frompos: pairs[k]=rara[k]
            else:
                pairs[':'.join([k, str(i)])]= rara[k]
                reader.fieldnames[i]=':'.join([k, str(i)])
        yield pairs


        #yield {':'.join([k, str(count)]): v for count, (k, v) in enumerate(rara.iteritems())}

def correct_fieldnames(reader):
    reader.fieldnames=[v if i<Columnnames.frompos else ':'.join([v, str(i)]) for i,v in enumerate(reader.fieldnames)]


def create_minimal_trio_bed(n,l):
    startpos = 0
    #print "#chrom\tchromStart\tchromEnd\tCode\tCall"
    for k, g  in itertools.groupby(l):
      ###!!!!!call len(list(g)) only once! every call probably next() method
      currentrange=len(list(g))
      print "%s\t%d\t%d\t%d\t%s"%(str(n),int(startpos),int(startpos+currentrange),int(k),names[k])
      startpos+=currentrange

def print_bed_header(): print "#chrom\tchromStart\tchromEnd\tCode\tCall"


with open(sys.argv[1], 'r') as csvfile:
  #chromosomes=ChromosomeCollection("1316_gDNA.GType","1316_1 PB2.GType",None)

  reader = csv.DictReader(csvfile,delimiter=';')
  metadata=MetaData(reader.fieldnames)
  #unique values for column header list
  correct_fieldnames(reader)


  #this will be replaced/extended by GUI/interactive menu
  gdna=metadata.get_gdna_by_sample("470")
  trios=metadata.get_trios_by_sample("470")
  ref=metadata.get_reference("470")
  #############################################
  chromosomes=ChromosomeCollection(gdna,ref,trios)
  countsNC=dict.fromkeys(reader.fieldnames[Columnnames.frompos:],0)
  countsAB=dict.fromkeys(reader.fieldnames[Columnnames.frompos:],0)
  alltogether=0


  for linecontent in reader:
      chrom=linecontent[Columnnames.chromosome]; name=linecontent[Columnnames.snpname]
      for k in countsAB.keys():
          chromosomes[chrom][k]=codes[linecontent[k]]
          if linecontent[k]=='AB':
              countsAB[k] += 1
          elif linecontent[k]=='NC':
              countsNC[k]+= 1

  records_num=reader.line_num-1
  #print "Processed %d records"  % (records_num)
  #print alltogether
  #print counts
  ####################################
  #Statistics -> show to the user
  ####################################
  '''
  print "%s\t%s\t%s" % ("Dataset","Calls","AB (proportion without NCs)")
  for k in sorted(countsAB.keys()):
    print "%s\t%.1f\t%.1f" % (k,round((1-countsNC[k]/float(records_num))*100,1),round(countsAB[k]/(float(records_num)-countsNC[k])*100,1))
  '''

  test=chromosomes['1'].phase()
  print_bed_header()
  chromosomes['1'].create_minimal_gdna_bed()
  create_minimal_trio_bed(test[0][0],test[0][1])
  create_minimal_trio_bed(test[1][0],test[1][1])
  create_minimal_trio_bed(test[2][0],test[2][1])
  ########################################
  #l1=chromosomes['1'].hetpos()
  #test=chromosomes['1'].phase()

  #l2=chromosomes['1']["1316_2 egg.GType"]


  #print len(l1)
  #print len(l2)
  #print len([item for item in l1 if item==1])


  #chromosomes['1'].dumptofile()


      #for k in (set(selection)& set(reader.fieldnames)):
      #    chromosomes[chrom][k]=codes[linecontent[k]]

  #print chromosomes
  #print (metadata._dataset)

  #print metadata._dataset

      #for rec in (set(selection)&set(linecontent.keys())): print linecontent[rec]

      #print linecontent
      #for rec in linecontent.keys()[4:]:
      #    print rec
  #sum(1 for row in reader if row["Grade"] == "A")
  #for row in csv.reader(input_file, delimiter=';'):
  #      grades[row[1]] += 1
  #print 'Number of A grades: %s' % grades['A']
  #print grades.most_common()

  #print counts_nc,linenr
  #print counts_ab,linenr


