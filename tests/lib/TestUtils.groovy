class TestUtils {

    def static floats_within_six_digits(Float float1, Float float2){
        def precision = 6;

        // Use an epsilon value based on the desired precision
        def epsilon = Math.pow(10, -precision)

        // Compare the absolute difference to the epsilon value
        return Math.abs(float1 - float2) < epsilon
    }

}
