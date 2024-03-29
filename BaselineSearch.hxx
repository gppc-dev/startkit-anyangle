#ifndef OPT_GPPC_BASELINE_SEARCH_HXX
#define OPT_GPPC_BASELINE_SEARCH_HXX

#include <vector>
#include <utility>
#include <limits>
#include <queue>
#include <stack>
#include <array>
#include <memory_resource>
#include <algorithm>
#include <cstdint>
#include <cassert>

namespace baseline
{

using std::uint32_t;
using std::size_t;

constexpr uint32_t COST_0 = 1'000;
constexpr uint32_t COST_1 = 1'414;
constexpr uint32_t COST_2 = 2'236;
using Point = std::pair<int,int>;
struct Node
{
	static constexpr uint32_t INV = std::numeric_limits<uint32_t>::max();
	static constexpr uint32_t FLOOD_FILL = INV-1;
	static constexpr uint32_t NO_PRED = INV-2;
	static constexpr uint32_t AMBIG_START = INV-3;
	static constexpr uint32_t AMBIG_NOGO = INV-4;
	// Node() noexcept : pred(INV), cost(INV)
	// { }
	uint32_t pred;
	uint32_t cost;
};
struct Grid
{
	size_t size() const noexcept { return cells.size(); }
	static Point pad(Point p) noexcept
	{
		return Point(p.first+1, p.second+1);
	}
	static Point unpad(Point p) noexcept
	{
		return Point(p.first-1, p.second-1);
	}
	uint32_t pack(Point p) const noexcept
	{
		assert(static_cast<uint32_t>(p.first) < width && static_cast<uint32_t>(p.second) < height);
		return static_cast<uint32_t>(p.second) * width + static_cast<uint32_t>(p.first);
	}
	Point unpack(uint32_t p) const noexcept
	{
		assert(width != 0 && p < size());
		return Point(static_cast<int>(p % width), static_cast<int>(p / width));
	}
	bool get(uint32_t p) const noexcept
	{
		return p < size() && cells[p];
	}
	bool get(Point p) const noexcept
	{
		return static_cast<uint32_t>(p.first) < width
		    && static_cast<uint32_t>(p.second) < height
			&& cells[pack(p)];
	}

	Grid(const std::vector<bool>& l_cells, int l_width, int l_height) :
		 width(static_cast<uint32_t>(l_width+2))
		,height(static_cast<uint32_t>(l_height+2))
	{
		// pad grid cells
		cells.assign(width * height, false);
		for (int y = 0; y < l_height; ++y)
		for (int x = 0; x < l_width; ++x) {
			cells[pack(Point(x+1, y+1))] = l_cells[y * l_width + x];
		}
		ambig.assign(cells.size(), false);
		uint32_t crd = 0;
		for (uint32_t y = 1; y < height; ++y) {
			assert(unpack(crd) == Point(0, y-1));
			std::uint8_t r1 = static_cast<std::uint8_t>(cells[crd]);
			std::uint8_t r2 = static_cast<std::uint8_t>(cells[crd+width]);
			for (uint32_t x = 1; x < width; ++x) {
				++crd;
				r1 = ((r1 << 1) | static_cast<std::uint8_t>(cells[crd])) & 0b11;
				r2 = ((r2 << 1) | static_cast<std::uint8_t>(cells[crd+width])) & 0b11;
				switch (r1 | (r2 << 2)) {
				case 0b1001:
				case 0b0110:
					ambig[crd+width] = true;
					break;
				default:
					break;
				}
			}
			++crd;
		}
	}

	uint32_t width;
	uint32_t height;
	std::vector<bool> cells;
	std::vector<bool> ambig;
	std::vector<Node> nodes;
};

void path_to_root(const Grid& grid, Point start, Point target, std::vector<Point>& out);
void setup_grid(Grid& grid);

struct SpanningTreeSearch : Grid
{
	SpanningTreeSearch(const std::vector<bool>& l_cells, int l_width, int l_height) : Grid(l_cells, l_width, l_height)
	{
		setup_grid(*this);
	}
	std::array<std::vector<Point>, 2> path_parts;
	const std::vector<Point>& get_path() const noexcept { return path_parts[0]; }
	// bool search found a path
	bool search(Point s, Point g)
	{
		s = pad(s);
		g = pad(g);
		std::array<uint32_t, 2> nodeid{{pack(s), pack(g)}};
		if (nodes[nodeid[0]].pred == Node::INV || nodes[nodeid[1]].pred == Node::INV)
			return false;
		if (nodeid[0] == nodeid[1]) {
			// zero path case
			path_parts[0].assign(2, unpad(s));
			return true;
		}
		path_parts[0].clear(); path_parts[1].clear();
		while (true) {
			int progressId = 0;
			if (auto c0 = nodes[nodeid[0]].cost, c1 = nodes[nodeid[1]].cost; c0 == c1) {
				// same dist, check if same root
				if (nodeid[0] == nodeid[1]) {
					path_parts[0].push_back(unpad(unpack(nodeid[progressId])));
					break; // found least common ancestor
				}
				if (c0 == 0)
					return false; // tree root's are different, no path
			} else if (c1 > c0) {
				progressId = 1; // nodeid[1] is longer thus process it first
			}
			path_parts[progressId].push_back(unpad(unpack(nodeid[progressId])));
			nodeid[progressId] = nodes[nodeid[progressId]].pred;
		}
		// form path
		path_parts[0].insert(path_parts[0].end(), path_parts[1].rbegin(), path_parts[1].rend());
		return true;
	}
};

void flood_fill(Grid& grid, std::pmr::vector<Point>& out, uint32_t origin, std::pmr::memory_resource* res)
{
	assert(origin < grid.nodes.size() && grid.nodes[origin].pred == Node::INV);
	out.clear();
	std::stack<Point, std::pmr::vector<Point>> Q(res);
	auto&& push_queue = [&grid,&Q] (Point p, int dx, int dy) {
		p.first += dx; p.second += dy;
		if (grid.get(p)) {
			if (auto& node = grid.nodes[grid.pack(p)]; node.pred == Node::INV) {
				node.pred = Node::FLOOD_FILL;
				Q.push(p);
			}
		}
	};
	grid.nodes[origin].pred = Node::FLOOD_FILL;
	Q.push(grid.unpack(origin));
	while (!Q.empty())
	{
		Point p = Q.top(); Q.pop();
		out.push_back(p);
		push_queue(p, 1, 0);
		push_queue(p, -1, 0);
		push_queue(p, 0, 1);
		push_queue(p, 0, -1);
	}
}

void dijkstra(Grid& grid, uint32_t origin, std::pmr::memory_resource* res)
{
	// first = dist, second = node-id
	using node_type = std::pair<uint32_t,uint32_t>;
	std::priority_queue<node_type, std::pmr::vector<node_type>, std::greater<node_type>> Q(res);
	auto try_push = [&grid,&Q](uint32_t node, int dx, int dy, uint32_t cost, bool push_ambig) {
		uint32_t newNode = static_cast<uint32_t>( static_cast<int>(node) + dy * grid.width + dx );
#ifndef NDEBUG
		Point p = grid.unpack(node);
		Point q = grid.unpack(newNode);
		assert(static_cast<uint32_t>(p.first + dx) < grid.width && static_cast<uint32_t>(p.second + dy) < grid.height);
		assert(p.first + dx == q.first && p.second + dy == q.second);
#endif
		bool ambig = grid.ambig[newNode];
		if ( (ambig & !push_ambig) )
			return;
		Node& N = grid.nodes[newNode];
		if (cost < N.cost) {
			N.pred = node;
			N.cost = cost;
			if (!ambig)
				Q.emplace(cost, newNode);
		}
	};
	Q.emplace(0, origin);
	grid.nodes[origin].cost = 0;
	grid.nodes[origin].pred = Node::NO_PRED;
	while (!Q.empty()) {
		auto [cost, node] = Q.top(); Q.pop();
		if (cost != grid.nodes[node].cost)
			continue; // skip
		Point p = grid.unpack(node);
		// push successors
		uint32_t mask = 0;
		for (int i = 0, dy = -2; dy < 2; dy++)
		for (int dx = -2; dx < 2; dx++) {
			mask |= static_cast<uint32_t>(grid.get( Point(p.first + dx, p.second + dy) )) << i++;
		}
		mask = ~mask; // 1 = non-trav, 0 = trav
		// 0123
		// 4567
		// 89AB
		// CDEF

		// ambig only first
		// E
		if ( (mask & 0b0000'0100'0100'0000u) != 0b0000'0100'0100'0000u ) try_push(node, 1, 0, cost + COST_0, false);
		// S
		if ( (mask & 0b0000'0110'0000'0000u) != 0b0000'0110'0000'0000u ) try_push(node, 0, 1, cost + COST_0, false);
		// SE
		if ( (mask & 0b0000'0100'0000'0000u) == 0 ) try_push(node, 1, 1, cost + COST_1, false);
		// S-SE
		if ( (mask & 0b0100'0100'0000'0000u) == 0 ) try_push(node, 1, 2, cost + COST_2, false);
		// E-SE
		if ( (mask & 0b0000'1100'0000'0000u) == 0 ) try_push(node, 2, 1, cost + COST_2, false);

		if ( (mask & 0b0000'0110'0110'0000u) != 0b0000'0010'0100'0000u ) {
			// not ambig point, permit all expansions
			// N
			if ( (mask & 0b0000'0000'0110'0000u) != 0b0000'0000'0110'0000u ) try_push(node, 0, -1, cost + COST_0, true);
			// W
			if ( (mask & 0b0000'0010'0010'0000u) != 0b0000'0010'0010'0000u ) try_push(node, -1, 0, cost + COST_0, true);
			// NE
			if ( (mask & 0b0000'0000'0100'0000u) == 0 ) try_push(node, 1, -1, cost + COST_1, false);
			// NW
			if ( (mask & 0b0000'0000'0010'0000u) == 0 ) try_push(node, -1, -1, cost + COST_1, true);
			// SW
			if ( (mask & 0b0000'0010'0000'0000u) == 0 ) try_push(node, -1, 1, cost + COST_1, false);
			// N-NE
			if ( (mask & 0b0000'0000'0100'0100u) == 0 ) try_push(node, 1, -2, cost + COST_2, false);
			// N-NW
			if ( (mask & 0b0000'0000'0010'0010u) == 0 ) try_push(node, -1, -2, cost + COST_2, true);
			// S-SW
			if ( (mask & 0b0010'0010'0000'0000u) == 0 ) try_push(node, -1, 2, cost + COST_2, false);
			// E-NE
			if ( (mask & 0b0000'0000'1100'0000u) == 0 ) try_push(node, 2, -1, cost + COST_2, false);
			// W-NW
			if ( (mask & 0b0000'0000'0011'0000u) == 0 ) try_push(node, -2, -1, cost + COST_2, true);
			// W-SW
			if ( (mask & 0b0000'0011'0000'0000u) == 0 ) try_push(node, -2, 1, cost + COST_2, false);
		}
	}
}

void setup_grid(Grid& grid)
{
	grid.nodes.assign(grid.size(), Node{Node::INV, Node::INV});
	std::pmr::unsynchronized_pool_resource vector_res;
	std::pmr::vector<Point> cluster;
	struct Dist {
		bool operator()(Point q, Point p) const noexcept {
			return dist(q, centre) < dist(p, centre);
		}
		static int dist(Point q, Point p) noexcept {
			return std::abs(q.first - p.first) + std::abs(q.second - p.second);
		}
		Point centre;
	};
	for (uint32_t i = 0, ie = grid.size(); i < ie; ++i) {
		if (grid.cells[i] && grid.nodes[i].pred == Node::INV) {
			// new cluster
			flood_fill(grid, cluster, i, &vector_res);
			assert(!cluster.empty());
			std::uint64_t sumx = 0, sumy = 0;
			for (Point p : cluster) {
				sumx += p.first; sumy += p.second;
			}
			Point cluster_centre(static_cast<int>(sumx / cluster.size()), static_cast<int>(sumy / cluster.size()));
			uint32_t cluster_id = grid.pack( *std::min_element(cluster.begin(), cluster.end(), Dist{cluster_centre}) );
			dijkstra(grid, cluster_id, &vector_res);
			assert(std::all_of(cluster.begin(), cluster.end(), [&grid] (Point q) { return grid.nodes.at(grid.pack(q)).pred != Node::FLOOD_FILL; }));
		}
	}
}

} // namespace baseline

#endif
