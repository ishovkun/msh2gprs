int discr_triplet_basic::
find_triplet(const Element & c,
             const Element & include,
             const Cell & primary,
             Storage::real v[3],
             Storage::real dir,
             Storage::real area,
             HandleType t[3],
             Storage::real coefs[3],
             Storage::real & in_rhs,
             Storage::real p_mult[3],
             const Element & other,
             bool positive_sum,
             bool print)
{
  const int kmax = 2048;
  int k = 0, best_wgt = INT_MAX, tri_wgt;

  Element elem[kmax];
  Storage::real rhs[kmax];
  Storage::real mult[kmax];
  Storage::real best_coef, max_coef;
  if (!primary.isValid())
  {
    std::cout << "Primary cell was not specified!" << std::endl;
    exit(-1);
  }

  std::vector<list_t> recs;
  int cid = c->LocalID();
  coord v0 = coord(v)*dir*area, unit_v0(v0.normalize());

  in_rhs = 0;
  t[0] = t[1] = t[2] = InvalidHandle();
  p_mult[0] = p_mult[1] = p_mult[2] = 0.0;
  coefs[0] = coefs[1] = coefs[2] = 0.0;

  if (length(v) < 1.0e-12) // no vector
    return 1;

  rMatrix mv = rMatrix::FromVector(v0.data(), 3), A(3, 3), coef(1, 3);
  rMatrix x0(3, 1), b(3, 1);
  c.Centroid(x0.data());

  start_search(primary);

  if (include->isValid())
  {
    rMatrix x1(3, 1), v(3, 1);
    if (include->GetElementType() == FACE)
      FindFaceCentroid(include->getAsFace(), x1.data());
    else
      include->Centroid(x1.data());
    v = x1 - x0;
    coord vr(v.data());
    mult[k] = 1;
    rhs[k] = 0;
    elem[k] = include;
    recs.push_back(list_t(vr, layer_wgt, k++));
  }

  rMatrix U, S, V;
  rMatrix Perm = rMatrix::FromTensor(primary.RealArray(K).data(), primary.RealArray(K).size());
  do {
    perform_search(primary, recs, x0, elem, mult, rhs, k, default_extend, print);

    if (!recs.empty())
    {
      bool necessery = include.isValid();
      if (layer_wgt > 2 && necessery) necessery = false;
      qsort(recs.begin() + (necessery ? 1 : 0), recs.end());

   reprint:
      int kk = -1;
      for (int q = 0; q < (int)recs.size() && kk == -1; ++q)
        if ((recs[q].p.normalize() - v0.normalize()).length() < 1.0e-8) kk = q;

      {
        if (layer_wgt > 4) break;
        bool test2d = true;

        {
          int l = -1;
          int a = -1;
          for (int k = 0; k < (int)recs.size() && a == -1; ++k) if ((recs[k].p).length() > 1.0e-6) a = k;
          for (int k = 0; k < (int)recs.size() && l == -1; ++k) if ((recs[a].p*recs[k].p).length() > 1.0e-6) l = k;
          coord nrm = recs[a].p*recs[l].p;
          nrm = nrm.normalize();

          if (fabs(nrm^v0) > 1.0e-7)
            test2d = false;

          for (int k = 0; k < (int)recs.size() && test2d; ++k)
            if (fabs(recs[k].p^nrm) > 1.0e-7)
            {
              if (print) printf("recs[%d].p %g %g %g recs[%d].p^nrm %g\n", k, recs[k].p[0], recs[k].p[1], recs[k].p[2], k, recs[k].p^nrm);
              test2d = false;
            }
        }

        Storage::real diag[3] = { 1, 1, 1 };
        if (test2d) diag[2] = 0;
        int dims = test2d ? 2 : 3, v[3];
        for (int d = 0; d < dims; ++d) v[d] = d;
        bool accepted, skip;
        do
        {
          tri_wgt = 0;
          for (int d = 0; d < dims; ++d)
            tri_wgt += recs[v[d]].wgt;

          if (tri_wgt <= best_wgt)
          {
            skip = false;
            for (int k = 0; k < dims && !skip; ++k)
            {
              for (int l = k + 1; l < dims && !skip; ++l)
                if (vec_colinear(recs[v[k]].p.data(), recs[v[l]].p.data())) skip = true;
            }
            if (!skip)
            {
              std::fill(A.data(), A.data() + 9, 0.0);
              double ln = 0;
              double cmults[3] = {0,0,0};
              for (int l = 0; l < dims; ++l) //go over nodes
              {
                int e = recs[v[l]].k;
                const coord & p = recs[v[l]].p;
                cmults[l] = mult[e];
                A(l, 0) = p[0];//*dist[e];
                A(l, 1) = p[1];//*dist[e];
                A(l, 2) = p[2];//*dist[e];
                ln += (p.normalize() - unit_v0).length()*p.length();
              }


              if( print )
              {
                std::cout << "Initial: " << std::endl;
                A.Print();
              }

              Threshold(A, 1.0e-12);

              if( print )
              {
                std::cout << "Threshold: " << std::endl;
                A.Print();
              }

              //std::pair<rMatrix, bool> inv = A.Invert();

              A.SVD(U, S, V);
              double Smax, Smin;
              Smax = Smin = S(0,0);
              for (int l = 0; l < dims; ++l)
              {
                if (S(l, l) > 0.0)
                {
                  Smin = S(l,l);
                  S(l, l) = 1.0 / S(l, l);
                }
                else S(l, l) = 0.0;
              }

              if( print )
              {
                std::cout << "S: " << std::endl;
                S.Print();
                std::cout << "U: " << std::endl;
                U.Print();
                std::cout << "V: " << std::endl;
                V.Print();
              }

              coef = mv.Transpose()*V*S*U.Transpose();

              bool nonnegative = true;
              if (positive_sum)
              {
                nonnegative &= (coef(0, 0)*cmults[0] + coef(0, 1)*cmults[1] + coef(0, 2)*cmults[2]) > 1.0e-15;
              }
              else
              {
                nonnegative &= coef(0, 0) > -1.0e-6;
                nonnegative &= coef(0, 1) > -1.0e-6;
                nonnegative &= coef(0, 2) > -1.0e-6;
                nonnegative &= (coef(0, 0)*cmults[0] + coef(0, 1)*cmults[1] + coef(0, 2)*cmults[2]) >= 0;// -1.0e-15;
              }

              if (nonnegative && (mv.Transpose()-coef*A).FrobeniusNorm() < 1.0e-4 )
              {

                {
                  max_coef = ln*Smax/Smin;
                  if (best_wgt == INT_MAX || max_coef < best_coef)
                  {
                    best_wgt = tri_wgt;
                    best_coef = max_coef;
                    in_rhs = 0.0;
                    double sum = 0, msum = 0;
                    for (int l = 0; l < dims; ++l)
                    {
                      coefs[l] = coef(0, l);
                      p_mult[l] = mult[recs[v[l]].k] * coefs[l];
                      sum += coefs[l];
                      msum += p_mult[l];
                      in_rhs += rhs[recs[v[l]].k] * coefs[l];
                      t[l] = elem[recs[v[l]].k].GetHandle();
                    }
                  }
                  // else if (print) printf("discarded\n");
                }
              // } else if( print ) printf("discarded, nonnegative or residual\n");
            } //else if (print) printf("system not solved\n");
          }
          // else if (print) printf("collinear\n");
        }
        accepted = false;
        for (int q = dims - 1; q >= 0 && !accepted; q--)
        {
          while (v[q] < (int)recs.size() && !accepted)
          {
            v[q]++;
            accepted = v[q] < (int)recs.size();
            for (int j = q + 1; j < dims; ++j)
            {
              v[j] = v[j - 1] + 1;
              if (v[j] >= (int)recs.size())
                accepted = false;
            }
          }
        }
        if (necessery && v[0] != 0) accepted = false;
      } while (accepted);
    }
  }
} while (best_wgt == INT_MAX && can_search());

  end_search();
  return best_wgt != INT_MAX;
}
