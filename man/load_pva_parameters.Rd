\name{loadParameters}
\alias{loadParameters}
\title{PVAInvasR Parameter Loader}
\usage{
loadParameters(csv_path)
}
\description{
Reads the .csv file located at `csv_path`, and put these into a form that can be accessed by PVAInvasR. Population and control parameters are stored in a list.

A template of the parameters needed by PVAInvasR can be downloaded to your computer using the function pvaTemplate().
}
\arguments{

\item{csv_path}{A string or file.path identifying the path to a filled-in parameter .csv.}

}
\value{
A list of loaded parameters in PVAInvasR foramt.
}
\examples{
pva_parameters <- loadParameters(csv_path)
}
