point_t
voxel_t::data_gradient(sdf_t * canon){
    return sdf->distance_gradient(p + u) * (sdf->distance(p + u) - canon->distance(p));
}

point_t
voxel_t::killing_gradient(float gamma){
    matrix_t jacobian = matrix_t::jacobian(sdf,  p + u);
    
    std::vector<float> jt_vec = jacobian.transpose().stack();
    std::vector<float> j_vec = jacobian.stack();
    std::vector<float> vec;
    for (int i = 0; i < 9; i++){
        vec.push_back(jt_vec[i] + gamma * j_vec[i]);
    }

    matrix_t h_u[3]; //TODO
   
    point_t result;
    for (int i = 0; i < 9; i++){
        result += point_t(
            h_u[i / 3].get(i % 3, 0) * vec[i]
            h_u[i / 3].get(i % 3, 1) * vec[i]
            h_u[i / 3].get(i % 3, 2) * vec[i]
        );
    }

    return 2 * result;
}

point_t
voxel_t::level_set_gradient(float epsilon){
    point_t g = sdf->distance_gradient(p + u);
    float scale = (g.length() - 1) / (g.length() + epsilon);
    
    matrix_t h = matrix_t::hessian(sdf, p + u);
                  
    return scale * (h * g);
}
