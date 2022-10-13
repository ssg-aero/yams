// #include <vector>
// #include <array>
// #include <algorithm>
// #include <execution>
// namespace yams
// {

//     template<typename T>
//     using vector = std::vector<T>;
//     template<typename T>
//     using point = std::array<T,2>;

//     template<typename T>
//     class bladeToBladeCurvatureSolver
//     {
//     private:
//         vector<T> X;
//         vector<T> Y;
//     public:
//         bladeToBladeCurvatureSolver(const vector<point<T>> &pts, size_t ni)
//         {
//             X.resize(pts.size());
//             Y.resize(pts.size());


//             std::transform(
//                 execution::par,
//                 pts.begin(), pts.end(),X.begin(),[](const auto &pt){return pt[0];}
//             )

//             std::transform(
//                 execution::par,
//                 pts.begin(), pts.end(),Y.begin(),[](const auto &pt){return pt[1];}
//             )
//         }

//     };

// }