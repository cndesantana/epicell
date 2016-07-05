#ifndef CELL_CLUSTER_IO_H
#define CELL_CLUSTER_IO_H

#include <fstream>
#include <string>
#include <vector>
#include "io.h"
#include "cellIo.h"
#include "core/epcDebug.h"
#include "division/cellDivisionHelpers.h"


namespace epc {

template<int d>
class VTKwriter {
public:
	VTKwriter(const CellCluster<d> &cellCluster_, std::string fname) : 
		cellCluster(cellCluster_), fout(fname.c_str()), firstCellField(true), firstVertexField(true)
	{ 
		numVerts = 0; // total number of vertices (redundant counting shared vertices are counted more than once...)
		for (auto c: cellCluster.getCells()) {
			numVerts += c->getEdges().size();
		}

		writeHeader();
		writeVertices();
		writePolygons();
	}

	~VTKwriter() {
		fout.close();
	}

	/// add a scalar value to each cell 
	void addScalarToCell(const std::string &sName, std::function<std::vector<double> (const CellCluster<d> &cellCluster)> &foo) {
		if (firstCellField) {
			fout << "CELL_DATA " << cellCluster.getCells().size() << std::endl;
			firstCellField = false;
		}
		fout << "SCALARS "+sName+" float" << std::endl;
		fout << "LOOKUP_TABLE default" << std::endl;
		auto vals = foo(cellCluster);
		for (auto v : vals) {
			fout << v << " ";
		}
		fout << std::endl;
	}

	/// add a scalar value to each vertex position
	void addScalarToVertex(const std::string &sName, std::function<std::vector<double> (const CellCluster<d> &cellCluster)> &foo) {
		if (firstVertexField) {
			fout << "POINT_DATA " << numVerts << std::endl;
			firstVertexField = false;
		}
		fout << "SCALARS "+sName+" float" << std::endl;
		fout << "LOOKUP_TABLE default" << std::endl;
		auto vals = foo(cellCluster);
		for (auto v : vals) {
			fout << v << " ";
		}
		fout << std::endl;
	}

	/// add a vector value to each cell 
	void addVectorToCell(const std::string &vName, std::function<std::vector<Array<double,d>> (const CellCluster<d> &cellCluster)> &foo) {
		if (firstCellField) {
			fout << "CELL_DATA " << cellCluster.getCells().size() << std::endl;
			firstCellField = false;
		}
		fout << "VECTORS "+vName+" float" << std::endl;
		auto vals = foo(cellCluster);
		for (auto v : vals) {
			fout << v << " 0 ";
		}
		fout << std::endl;
	}

	/// add a vector value to each vertex position 
	void addVectorToVertex(const std::string &vName, std::function<std::vector<Array<double,d>> (const CellCluster<d> &cellCluster)> &foo) {
		if (firstVertexField) {
			fout << "POINT_DATA " << numVerts << std::endl;
			firstVertexField = false;
		}
		fout << "VECTORS "+vName+" float" << std::endl;
		auto vals = foo(cellCluster);
		for (auto v : vals) {
			fout << v << " 0 ";
		}
		fout << std::endl;
	}

private:
	/// mandatory header for the vtk file
	void writeHeader() {
		fout << "# vtk DataFile Version 3.0" << std::endl;
		fout << "vtk output" << std::endl;
		fout << "ASCII" << std::endl;
		fout << "DATASET POLYDATA" << std::endl;
	}

	/// write the position of each vertex (again redundant positioning, same vertex is written more than once)
	void writeVertices() {
		fout << "POINTS " << numVerts << " float" << std::endl;
	    
		for (auto c: cellCluster.getCells()) {
			const auto verts = vertexHelpers<d>::sortVerticesCounterClockwise(cellHelpers<d>::getVertices(*c));
			for (auto v: verts) {
				fout << v->getPosition() << " 0" << std::endl;
			}
		}
	}

	/// writes the connectivity of the polygons 
	void writePolygons() {
		fout << "POLYGONS " << cellCluster.getCells().size() << " " << cellCluster.getCells().size()+numVerts << std::endl;
		unsigned long offset = 0;
		for (auto c: cellCluster.getCells()) {
			unsigned size = c->getEdges().size();
			fout << size;
			for (unsigned iA = 0; iA < size; ++iA) {
				fout << " " << offset+iA;
			}
			offset += size;
			fout << std::endl;
		}
	}

private:
	CellCluster<d> cellCluster;
	std::ofstream fout;
	bool firstCellField,firstVertexField;
	unsigned long numVerts;
};


template<int d>
struct cellClusterIO {

static void writeEdges(const CellCluster<d> &cellCluster, std::ofstream &fout) {
//	unsigned long idCell = 0;
	unsigned long nVerts = 0;
	for (auto c: cellCluster.getCells()) {
		auto idCell = c->getId();
		auto listOfVerts = cellHelpers<d>::getVertices(*c);
		const auto verts = vertexHelpers<d>::sortVerticesCounterClockwise(cellHelpers<d>::getVertices(*c));
		nVerts = verts.size();
	    	for (unsigned iV = 0; iV < verts.size(); ++iV) {
//			CellDivisionHelpers<d>::getEdgesByVertices(cellCluster.getEdges(),verts[iV],verts[(iV+1)%nVerts])->getId();	
			fout << idCell << " " << c->getGamma() << " " << verts[iV]->getPosition() << " " << verts[(iV+1)%nVerts]->getPosition() << std::endl;
		}
	}
}

static void writeRegularness(std::ofstream &fout,int iT,double regularness, double gammaNorm, double lambdaNorm, double sdgamma, double sdlambda){
	fout << iT << " " << gammaNorm << " " << lambdaNorm << " " <<  sdgamma << " " << sdlambda << " " << regularness << std::endl; 
}

static void writeDistributionOfSides(std::ofstream &fout,int iT, std::vector<double> vecOfProportions){
	fout << iT << " "; 
	for (auto iv: vecOfProportions){
            fout << iv << " ";
        } 
	fout << std::endl; 
}

static void writeConvergence(std::ofstream &fout, int iT, double gammaNorm, double lambdaNorm, double gamma, double lambda, int convergence) {

	fout << gammaNorm << " " << lambdaNorm << " " << gamma << " " << lambda << " " << convergence << std::endl; 

}

static void writeStatistics(const CellCluster<d> &cellCluster, std::ofstream &fout, int iT, double gammaNorm, double lambdaNorm, double sdgamma, double sdlambda) {
	auto nCells = cellCluster.getCells().size();
	for (auto c: cellCluster.getCells()) {
		auto idCell = c->getId();
		auto areaCell = c->computeArea();
		auto A0Cell = c->getA0();
		auto perimeterCell = c->computePerimeter();
		auto shapeCell = c->getEdges().size();

		fout << iT << " " << nCells << " " << idCell << " " << gammaNorm << " " << lambdaNorm << " " << sdgamma << " " << sdlambda << " " << areaCell << " " << A0Cell << " " << perimeterCell << " " << shapeCell << std::endl;
	}
}

static void writeDIST(std::string fname, long unsigned int iT, std::vector<double> vecOfProportions){
    std::ofstream fout(fname.c_str(),std::fstream::in | std::fstream::out | std::fstream::app);
    writeDistributionOfSides(fout,iT,vecOfProportions);
    fout.close();
}

static void writeREG(std::string fname, long unsigned int iT, double regularness, double gammaNorm, double lambdaNorm, double sdgamma, double sdlambda){
    std::ofstream fout(fname.c_str(),std::fstream::in | std::fstream::out | std::fstream::app);
    writeRegularness(fout,iT,regularness,gammaNorm,lambdaNorm,sdgamma,sdlambda);
    fout.close();
}

static void writeCVG(std::string fname, int iT, double gammaNorm, double lambdaNorm, double gamma, double lambda, int convergence){
    std::ofstream fout(fname.c_str(),std::fstream::in | std::fstream::out | std::fstream::app);
    writeConvergence(fout,iT,gammaNorm,lambdaNorm,gamma,lambda,convergence);
    fout.close();
}

/// writes cellCluster as a list of edges (Edges List Format) 
static void writeSTA(const CellCluster<d> &cellCluster, std::string fname, int iT, double gammaNorm, double lambdaNorm, double sdgamma, double sdlambda){
    std::ofstream fout(fname.c_str(),std::fstream::in | std::fstream::out | std::fstream::app);
    writeStatistics(cellCluster,fout,iT,gammaNorm,lambdaNorm,sdgamma,sdlambda);
    fout.close();
}

/// writes cellCluster as a list of edges (Edges List Format) 
static void writeELF(const CellCluster<d> &cellCluster, std::string fname){
    std::ofstream fout(fname.c_str(),std::fstream::out | std::fstream::trunc);
    writeEdges(cellCluster,fout);
    fout.close();
}
};// end cellClusterIO

} // end namespace epc

#endif
