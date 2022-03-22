FROM alpine:3.14 as gurobi
ARG GRB_VERSION=9.5.1
ARG GRB_SHORT_VERSION=9.6

WORKDIR /opt

RUN apk add --no-cache wget tar \
    && wget -v https://packages.gurobi.com/${GRB_SHORT_VERSION}/gurobi${GRB_VERSION}_linux64.tar.gz \
    && tar -xvf gurobi${GRB_VERSION}_linux64.tar.gz  \
    && rm -f gurobi${GRB_VERSION}_linux64.tar.gz \
    && mv -f gurobi* gurobi \
    && rm -rf gurobi/linux64/docs


FROM gradle:jdk17 as gradle

COPY --chown=gradle:gradle sources.tar.gz /opt/sources.tar.gz
WORKDIR /opt
RUN tar -xf sources.tar.gz \
   && cd /opt/java \
   && gradle distTar --no-daemon


FROM openjdk:17-slim

COPY --from=gurobi /opt/gurobi /opt/gurobi
COPY --from=gradle /opt/java/build/distributions/ssp-cvar-1.0.tar /app/
COPY --from=gradle /opt/LICENCE /app/LICENCE
COPY --from=gradle /opt/models /app/models
COPY --from=gradle /opt/python /app/
WORKDIR /app
RUN tar -xvf ssp-cvar-1.0.tar --strip-components=1 \
    && rm ssp-cvar-1.0.tar

ENV GUROBI_HOME="/opt/gurobi/linux64"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$GUROBI_HOME/lib"
ENV PATH="$PATH:$GUROBI_HOME/bin"