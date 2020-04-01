function xi = groupMult(xi1,xi2)
    % groupMult - SE3 lie Group multiplication
    %
    % xi1, xi2: input Lie algebra elements

    xi = transform_SE3_se3_homCoords(se3_SE3_homCoords(xi1)*se3_SE3_homCoords(xi2));
end

