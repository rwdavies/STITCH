// not ready to properly junk yet


//  //' @export
//  // [[Rcpp::export]]
// Rcpp::List rcpp_sample_multiple_paths(const int n_starts, const int n_its, const Rcpp::List& sampleReads, const int nReads, const arma::mat& eHaps_t, const double maxDifferenceBetweenReads, const int Jmax, const arma::vec pi, const arma::mat& transMatRate_t, const arma::mat& alphaMat_t, const arma::ivec& srp, const arma::ivec& sum_dosage_vec) {
//     //
//     // new variables go here
//     //
//     double n_save_iterations;
//     const int T = eHaps_t.n_cols;
//     int i_start, iRead, iHap, it, t;
//     double s1, s2, p1, p2, pHap1;    
//     arma::imat path(2, T);
//     arma::mat read_labels = arma::zeros(2, nReads);
//     arma::mat dosages = arma::zeros(n_starts, T);
//     //
//     // output variables go here
//     //
//     arma::mat eMatHapPH_t;
//     arma::vec pRgivenH1;
//     arma::vec pRgivenH2;
//     arma::mat eMatHap_t = rcpp_make_eMatHap_t(
//         sampleReads,
//         nReads,
//         eHaps_t,
//         maxDifferenceBetweenReads,
//         Jmax,
//         eMatHapPH_t,
//         pRgivenH1,
//         pRgivenH2
//     );
//     //
//     // initialize random numbers
//     //
//     double rand_uniform = 0;
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_real_distribution<> dis(0, 1);
//     //
//     // loop over starts
//     //
//     for(i_start = 0; i_start < n_starts; i_start++) {
//         // start with random labels
//         for(iRead=0; iRead < nReads; iRead++) {
//             rand_uniform = dis(gen);
//             if (rand_uniform > 0.5) {
//                 read_labels(0, iRead) = 1;
//             } else {
//                 read_labels(1, iRead) = 1;
//             }
//         }
//         for(it = 0; it < n_its; it++) {
//             // sample path
//             for(iHap = 0; iHap < 2; iHap++) {
//                 path.row(iHap) = rcpp_sample_path(read_labels.row(iHap), eMatHap_t, sampleReads, nReads, eHaps_t, maxDifferenceBetweenReads, Jmax, pi, transMatRate_t, alphaMat_t);
//             }
//             // sample labels here
//             for(iRead = 0; iRead < nReads; iRead ++) {
//                 s1 = path(0, srp(iRead));
//                 s2 = path(1, srp(iRead));
//                 p1 = 0.5 * eMatHap_t(s1, iRead);
//                 p2 = 0.5 * eMatHap_t(s2, iRead);
//                 pHap1 = p1 / (p1 + p2);
//                 rand_uniform = dis(gen);
//                 if (rand_uniform < pHap1) {
//                     read_labels(0, iRead) = 1;
//                     read_labels(1, iRead) = 0;
//                 } else {
//                     read_labels(0, iRead) = 0;
//                     read_labels(1, iRead) = 1;
//                 }
//             }
//             // can save dosages here
//             if (sum_dosage_vec(it) == 1) {
//                 for(iHap = 0; iHap < 2; iHap++) {
//                     for(t=0; t < T; t++) {
//                         dosages(i_start, t) = dosages(i_start, t) +      \
//                             eHaps_t(path(iHap, t), t);
//                     }
//                 }
//             }

//         } // end of iterations
//         //
//         // normalize dosage here
//         //
//         n_save_iterations = 0;
//         for(int i=0; i < n_its; i++)
//             n_save_iterations = n_save_iterations + sum_dosage_vec(i);
//         for(t=0; t < T; t++) {
//             dosages(i_start, t) = dosages(i_start, t) / n_save_iterations;
//         }
//     } // end of starts
//     return(wrap(Rcpp::List::create(
//                                    Rcpp::Named("dosages") = dosages,
//                                    Rcpp::Named("path") = path,
//                                    Rcpp::Named("read_labels") = read_labels
//                                    )));
// }
