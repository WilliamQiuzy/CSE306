// TD9

std::vector<Vector> embedding(std::vector<Vector> vertices, const std::vector<std::vector<int>> &adjacencyList, const std::vector<int> &boundaryVertices, int iterations) 
{
    int boundarySize = boundaryVertices.size();
    double perimeterLength = 0;

    for (int i = 0; i < boundarySize; i++) {
        perimeterLength += std::sqrt((vertices[boundaryVertices[i]] - vertices[boundaryVertices[(i + 1) % boundarySize]]).norm());
    }

    double cumulativeLength = 0;
    for (int i = 0; i < boundarySize; i++) {
        double angle = 2 * M_PI * cumulativeLength / perimeterLength;
        vertices[boundaryVertices[i]] = Vector(std::cos(angle), std::sin(angle), 0);
        cumulativeLength += std::sqrt((vertices[boundaryVertices[i]] - vertices[boundaryVertices[(i + 1) % boundarySize]]).norm());
    }

    for (int iteration = 0; iteration < iterations; iteration++) {
        std::vector<Vector> newVertices = vertices;

        for (size_t i = 0; i < vertices.size(); i++) {
            newVertices[i] = Vector(0, 0, 0);
            size_t neighborsCount = adjacencyList[i].size();

            for (int neighbor : adjacencyList[i]) {
                newVertices[i] = newVertices[i] + vertices[neighbor];
            }
            newVertices[i] = newVertices[i] / neighborsCount;
        }



        vertices = newVertices;
    }

    return vertices;
}