/*
Copyright (c) 2023 Grid-based Path Planning Competition and Contributors <https://gppc.search-conference.org/>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef GPPC_AA_VALIDATEPATH_HPP_INCLUDED
#define GPPC_AA_VALIDATEPATH_HPP_INCLUDED

#include "BresenhamRay.hpp"

namespace inx {

std::unique_ptr<BresenhamRay>& ValidatePath_data()
{
    static std::unique_ptr<BresenhamRay> rayShooter;
    return rayShooter;
}

// returns -1 if valid path, otherwise id of segment where invalidness was detetcted
// map is loaded only ONCE, if has changed, will error
// can call ValidatePath_data().reset() to allow new map
template <typename T>
int ValidatePath(const std::vector<bool>& map, int width, int height, const T& path)
{
    auto& rayShooter = ValidatePath_data();
    if (rayShooter == nullptr) {
        rayShooter = std::make_unique<BresenhamRay>();
        rayShooter->setGrid<true>(static_cast<size_t>(width), static_cast<size_t>(height), map);
    }
    return rayShooter->validPath(path);
}

} // namespace inx

#endif // GPPC_AA_VALIDATEPATH_HPP_INCLUDED
