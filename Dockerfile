FROM python:3.7
RUN apt-get update \
        && apt-get install -y python3-pip \
        # dependencies for jSVD
        && pip install numpy==1.16\
        && pip install scipy==1.6.0\
        && pip install autograd==1.3\
        && pip install pymanopt==0.2.5\
