% Calculate inertial parameters regressor of potential energy for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh2m2DE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:44
% EndTime: 2020-05-03 01:06:44
% DurationCPUTime: 0.12s
% Computational Cost: add. (58->28), mult. (91->30), div. (0->0), fcn. (77->8), ass. (0->19)
t52 = cos(qJ(1));
t60 = g(1) * t52;
t48 = sin(qJ(1));
t59 = g(2) * t48;
t46 = sin(qJ(3));
t47 = sin(qJ(2));
t56 = t47 * pkin(4);
t58 = g(3) * (t46 * pkin(5) + t56);
t57 = g(3) * t47;
t51 = cos(qJ(2));
t55 = t51 * pkin(4) + pkin(1);
t54 = pkin(2) + t55;
t50 = cos(qJ(3));
t53 = t50 * pkin(5) + t54;
t40 = t59 + t60;
t49 = cos(qJ(4));
t45 = sin(qJ(4));
t39 = -g(1) * t48 + g(2) * t52;
t1 = [0, 0, 0, 0, 0, 0, -t40, -t39, -g(3), 0, 0, 0, 0, 0, 0, 0, -t40 * t51 - t57, -t51 * g(3) + t40 * t47, t39, -pkin(1) * t40, 0, 0, 0, 0, 0, 0, -t40, -g(3), t39, -t55 * t60 - pkin(1) * t59 + (-t51 * t59 - t57) * pkin(4), 0, 0, 0, 0, 0, 0, -g(3) * t46 - t40 * t50, -t50 * g(3) + t40 * t46, t39, -g(3) * t56 - t40 * t54, 0, 0, 0, 0, 0, 0, -t40, -g(3), t39, -t40 * t53 - t58, 0, 0, 0, 0, 0, 0, -t45 * t39 - t49 * t40, -t39 * t49 + t45 * t40, -g(3), -t58 - t40 * (pkin(3) + t53);];
U_reg = t1;
