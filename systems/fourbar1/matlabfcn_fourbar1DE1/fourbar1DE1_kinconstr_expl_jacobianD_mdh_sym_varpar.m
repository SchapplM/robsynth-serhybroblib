% Jacobian time derivative of explicit kinematic constraints of
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% WD [5x1]
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)
% [ParkChoPlo1999] Park, FC and Choi, Jihyeon and Ploen, SR: Symbolic formulation of closed chain dynamics in independent coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function WD = fourbar1DE1_kinconstr_expl_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_kinconstr_expl_jacobianD_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_kinconstr_expl_jacobianD_mdh_sym_varpar: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_kinconstr_expl_jacobianD_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:54:31
% EndTime: 2020-04-24 19:54:31
% DurationCPUTime: 0.09s
% Computational Cost: add. (330->30), mult. (449->57), div. (6->4), fcn. (116->4), ass. (0->33)
t71 = pkin(1) ^ 2;
t67 = cos(qJ(1));
t85 = pkin(2) * t67;
t79 = -0.2e1 * pkin(1) * t85 + t71;
t89 = -pkin(3) - pkin(4);
t55 = (pkin(2) - t89) * (pkin(2) + t89) + t79;
t88 = -pkin(3) + pkin(4);
t56 = (pkin(2) - t88) * (pkin(2) + t88) + t79;
t81 = t55 * t56;
t72 = sqrt(-t81);
t53 = 0.1e1 / t72;
t77 = qJD(1) * t53;
t90 = 0.2e1 * pkin(2);
t66 = sin(qJ(1));
t86 = pkin(2) * t66;
t84 = (-t55 - t56) * pkin(1) * t86 * t77;
t70 = pkin(2) ^ 2;
t87 = pkin(1) * t70;
t83 = 0.1e1 / t81 * t84;
t61 = t70 + t79;
t59 = 0.1e1 / t61;
t82 = t53 * t59;
t80 = t66 * t72;
t78 = pkin(3) ^ 2 - pkin(4) ^ 2;
t76 = t71 * t90;
t75 = t59 * t83;
t74 = 0.1e1 / t61 ^ 2 * t66 * t77;
t65 = t66 ^ 2;
t63 = -pkin(1) * t67 + pkin(2);
t62 = pkin(1) - t85;
t58 = t61 - t78;
t57 = t61 + t78;
t1 = [0; -(t62 * t84 + (0.2e1 * t65 * t87 + (t58 * t67 + t80) * pkin(2)) * qJD(1)) * pkin(1) * t82 + (-pkin(1) * t75 + t74 * t76) * (t58 * t86 + t62 * t72); pkin(2) * (t63 * t84 + (t65 * t76 + (t57 * t67 + t80) * pkin(1)) * qJD(1)) * t82 + (pkin(2) * t75 - 0.2e1 * t74 * t87) * (pkin(1) * t66 * t57 + t63 * t72); (t66 * t83 + t67 * t77) * pkin(1) * t90; 0;];
WD = t1;
