% Calculate minimal parameter regressor of joint inertia matrix for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% MM_reg [((1+1)*1/2)x9]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:21
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = fourbar1TE_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_inertiaJ_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_inertiaJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:21:02
% EndTime: 2020-06-26 17:21:03
% DurationCPUTime: 0.14s
% Computational Cost: add. (190->55), mult. (312->81), div. (10->3), fcn. (48->4), ass. (0->37)
t35 = pkin(1) ^ 2;
t25 = cos(qJ(1));
t55 = pkin(2) * t25;
t46 = pkin(1) * t55;
t49 = -0.2e1 * t46 + t35;
t27 = pkin(3) - pkin(4);
t52 = (pkin(2) + t27) * (pkin(2) - t27);
t26 = pkin(3) + pkin(4);
t53 = (pkin(2) + t26) * (pkin(2) - t26);
t57 = (t49 + t53) * (t49 + t52);
t24 = sin(qJ(1));
t56 = pkin(1) * t24;
t33 = pkin(2) ^ 2;
t11 = t33 + t49;
t54 = 0.1e1 / t11 ^ 2 * t33;
t12 = pkin(1) * t25 - pkin(2);
t36 = sqrt(-t57);
t3 = t12 * t36;
t51 = t27 ^ 2 * t26 ^ 2;
t19 = t25 ^ 2;
t50 = t35 * t19;
t29 = pkin(4) ^ 2;
t48 = -t29 / 0.2e1 + t33;
t30 = pkin(3) ^ 2;
t47 = t29 - t30;
t45 = t30 / 0.2e1 + t48;
t44 = t33 - t47;
t43 = t54 / t57;
t42 = 0.1e1 / pkin(3) / t36 * t54;
t37 = pkin(1) * t35;
t34 = t35 ^ 2;
t32 = t33 ^ 2;
t28 = 0.2e1 * t33;
t14 = t29 + t30;
t2 = t3 + (t11 + t47) * t56;
t1 = t3 - (t44 + t49) * t56;
t4 = [1, 0, 0, -t2 ^ 2 * t43, -((0.4e1 * t45 * t50 - 0.4e1 * (t35 + t45) * t46 + t34 + (t28 + t47) * t35 + t33 * t44) * t36 + 0.2e1 * (t11 * t14 - t51) * t12 * t56) * t42, -0.2e1 * (0.6e1 * ((t33 - t30 / 0.6e1 - t29 / 0.6e1) * t35 + t32 + (-0.5e1 / 0.6e1 * t30 - 0.5e1 / 0.6e1 * t29) * t33 + t51 / 0.6e1) * t50 + 0.3e1 / 0.2e1 * t34 * t33 + (0.3e1 / 0.2e1 * t32 - t14 * t33 - t51 / 0.2e1) * t35 + t33 * t52 * t53 / 0.2e1 + (t24 * t27 * t26 * t3 - 0.3e1 * (t34 + (t28 - 0.2e1 / 0.3e1 * t30 - 0.2e1 / 0.3e1 * t29) * t35 + t32 + (-0.4e1 / 0.3e1 * t30 - 0.4e1 / 0.3e1 * t29) * t33 + t51 / 0.3e1) * t55) * pkin(1) + (-0.4e1 * (-t30 / 0.2e1 + t48) * t19 * t55 + t37 / 0.2e1) * t37) * t42, -t1 ^ 2 * t43, 0, 0;];
MM_reg = t4;
