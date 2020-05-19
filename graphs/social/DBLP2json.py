#!/usr/bin/python
# wget -N http://dblp.uni-trier.de/xml/dblp.xml.gz
# then run this script

import codecs, collections, json, gzip, os, sys, xml.sax
import io
import pickle
xml_filename = 'dblp.xml.gz'
json_gz_filename = 'dblp.json.gz'
tmp_filename = 'tmp_dblp.json.gz'
report_frequency = 10000

authorsEvents = {}
authorsAreas = {}

theory = ['stoc', 'soda', 'focs', 'icalp']
dmw = ['kdd', 'icdm', 'wsdm', 'www']
ivc = ['cvpr', 'iccv', 'eccv']
ml = ['nips', 'colt', 'icml']
idm = ['sigmod', 'pods', 'vldb']
net = ['sigcomm', 'infocom', 'mobicom', 'imc']

eventIndex = {}
for i, category in enumerate([theory, dmw, ivc, ml, idm, net]):
    eventIndex.update({event: i for event in category})


class DBLPHandler (xml.sax.ContentHandler):
    papertypes = set(['inproceedings', 'proceedings'])  # 'set (['article', 'book', 'inproceedings', 'incollection', 'www', 'proceedings', 'phdthesis', 'mastersthesis'])
    conferences = ['stoc', 'soda', 'focs', 'icalp', 'kdd', 'icdm', 'wsdm', 'www', 'cvpr', 'iccv' 'eccv', 'nips', 'colt', 'icml', 'sigmod', 'pods', 'vldb', 'sigcomm', 'infocom', 'mobicom', 'imc']

    def __init__(self, out):
        self.out = out
        self.paper = None
        self.authors = []
        self.year = None
        self.text = ''
        self.papercount = 0
        self.edgecount = 0

    def startElement(self, name, attrs):
        if name in self.papertypes:
            self.paper = str(attrs['key'])
            self.authors = []
            self.year = None
        elif name in ['author', 'year']:
            self.text = ''

    def endElement(self, name):
        if name == 'author':
            self.authors.append(self.text)
        if name == 'year':
            self.year = int(self.text.strip())
        elif name in self.papertypes:
            self.write_paper()
            self.paper = None

    def write_paper(self):
        if any(conference in self.paper for conference in self.conferences) and self.year >= 1998:
            self.addEventToAuthors()
            if self.papercount:
                self.out.write(',\n')
            self.papercount += 1
            self.edgecount += len(self.authors)
            json.dump([self.paper, self.authors, self.year], self.out)
            if self.papercount % report_frequency == 0:
                print('... processed %d papers, %d edges so far ...' % (self.papercount, self.edgecount))
                sys.stdout.flush()

    def characters(self, chars):
        self.text += chars

    def addEventToAuthors(self):
        for event in self.conferences:
            if event in self.paper:
                for author in self.authors:
                    if author in authorsEvents:
                        if event not in authorsEvents[author]:
                            authorsEvents[author].append(event)
                            authorsAreas[author][eventIndex[event]] += 1
                    else:
                        authorsEvents[author] = [event]
                        authorsAreas[author] = [0, 0, 0, 0, 0, 0]
                        authorsAreas[author][eventIndex[event]] += 1


def force():
    print('** Parsing XML...')
    f1 = io.open('test', 'w')
    f1.close()
    xmlfile = gzip.GzipFile(xml_filename, 'r')
    out = gzip.GzipFile(tmp_filename, 'w')
    out.write('[\n')
    dblp = DBLPHandler(out)
    parser = xml.sax.parse(xmlfile, dblp)
    out.write('\n]\n')
    out.close()
    os.rename(tmp_filename, json_gz_filename)
    with open('authorsEvents.pickle', 'wb') as f1:
        pickle.dump(authorsEvents, f1, protocol=pickle.HIGHEST_PROTOCOL)
    with open('authorsAreas.pickle', 'wb') as f2:
        pickle.dump(authorsAreas, f2, protocol=pickle.HIGHEST_PROTOCOL)

    print('-- %d papers, %d edges' % (dblp.papercount, dblp.edgecount))


def main(parse_args=False):
    try:
        need = (os.stat(xml_filename).st_mtime >= os.stat(json_gz_filename).st_mtime)
    except OSError:
        need = True
    if parse_args and len(sys.argv) > 1:
        need = True
    need = True
    if need:
        force()


def bropen():
    main()
    return gzip.GzipFile(json_gz_filename, 'r')


def papers():
    for line in bropen():
        if line.strip() in '[]': continue
        line = line.rstrip().rstrip(',')
        yield json.loads(line)


if __name__ == '__main__':
    main(True)
