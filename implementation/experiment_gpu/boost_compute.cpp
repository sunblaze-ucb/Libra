#include <vector>
#include <algorithm>

#include <boost/compute/algorithm/transform.hpp>
#include <boost/compute/container/vector.hpp>
#include <boost/compute/functional/math.hpp>

namespace compute = boost::compute;

int main()
{
    // get default device and setup context
    compute::device device = compute::system::default_device();
    compute::context context(device);
    compute::command_queue queue(context, device);

    // generate random data on the host
    std::vector<float> host_vector(10000);
    std::generate(host_vector.begin(), host_vector.end(), rand);

    // create a vector on the device
    compute::vector<float> device_vector(host_vector.size(), context);

    // transfer data from the host to the device
    compute::copy(
        host_vector.begin(), host_vector.end(), device_vector.begin(), queue
    );

    // calculate the square-root of each element in-place
    compute::transform(
        device_vector.begin(),
        device_vector.end(),
        device_vector.begin(),
        compute::sqrt<float>(),
        queue
    );

    // copy values back to the host
    compute::copy(
        device_vector.begin(), device_vector.end(), host_vector.begin(), queue
    );

    return 0;
}
