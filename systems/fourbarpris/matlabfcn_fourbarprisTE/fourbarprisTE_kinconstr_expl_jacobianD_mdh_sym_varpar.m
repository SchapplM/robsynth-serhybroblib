% Jacobian time derivative of explicit kinematic constraints of
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% WD [5x1]
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)
% [ParkChoPlo1999] Park, FC and Choi, Jihyeon and Ploen, SR: Symbolic formulation of closed chain dynamics in independent coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function WD = fourbarprisTE_kinconstr_expl_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_kinconstr_expl_jacobianD_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_kinconstr_expl_jacobianD_mdh_sym_varpar: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_kinconstr_expl_jacobianD_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 08:59:52
% EndTime: 2020-05-07 08:59:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (251->21), mult. (118->26), div. (6->3), fcn. (16->2), ass. (0->22)
t50 = qJ(1) + pkin(3);
t62 = -pkin(2) - t50;
t44 = pkin(1) - t62;
t61 = -pkin(2) + t50;
t45 = pkin(1) - t61;
t46 = pkin(1) + t61;
t47 = pkin(1) + t62;
t65 = t46 * t47;
t60 = t45 * t65;
t59 = t44 * t60;
t39 = ((-t59) ^ (-0.1e1 / 0.2e1));
t68 = (qJD(1) * t39);
t67 = (-t60 + (t65 + (t46 - t47) * t45) * t44) / t59 * t68;
t48 = 0.1e1 / t50;
t66 = (t50 * qJD(1) * t48);
t64 = (pkin(1) ^ 2 - pkin(2) ^ 2);
t63 = (qJD(1) / t50 ^ 2);
t58 = t48 * t67 / 0.2e1;
t57 = -qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t42 = t57 + t64;
t41 = -t57 + t64;
t1 = [t42 * t58 + ((-t42 * t63 - 2 * t66) * t39); 0; -t50 * t67 - (2 * t68); t41 * t58 + ((-t41 * t63 + 2 * t66) * t39); 0;];
WD = t1;
