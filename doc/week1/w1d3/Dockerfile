# build me with `docker build -t naturalis/mebioda .`
# run me with `docker run -it naturalis/mebioda`
# publish me with `docker login` and `docker push naturalis/mebioda`
FROM ubuntu:18.04
RUN apt-get update && apt-get install -y python-pip curl muscle mafft nano
RUN pip install biopython simplejson
RUN curl -L -o /usr/local/src/raxml.tar.gz https://github.com/stamatak/standard-RAxML/archive/v8.2.12.tar.gz
RUN cd /usr/local/src && tar xvzf raxml.tar.gz
RUN cd /usr/local/src/standard-RAxML-8.2.12 && make -f Makefile.gcc
ENV PATH "${PATH}:/usr/local/src/standard-RAxML-8.2.12"
CMD [ "/bin/bash" ]

