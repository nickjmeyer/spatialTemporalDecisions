#include "tuneParam.hpp"

bool TuneParam::getEdgeToEdge() const {
	return this->edgeToEdge;
}

void TuneParam::setEdgeToEdge(const bool edgeToEdge) {
	this->edgeToEdge = edgeToEdge;
}
