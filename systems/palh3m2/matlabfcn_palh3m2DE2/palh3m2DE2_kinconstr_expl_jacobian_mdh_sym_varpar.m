% Jacobian of explicit kinematic constraints of
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% W [12x4]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = palh3m2DE2_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:12:57
% EndTime: 2020-05-07 02:12:57
% DurationCPUTime: 0.16s
% Computational Cost: add. (754->0), mult. (1220->0), div. (53->0), fcn. (1828->0), ass. (0->1)
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, -1, -1, 0; 0, 0, 0, 1; 0, 1, 0, 0; 0, -1, 0, 0; 0, 1, 1, 0; 0, 1, 0, 0; 0, 1, 1, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
W = t1;
