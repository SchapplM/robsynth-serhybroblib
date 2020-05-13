% Calculate inertial parameters regressor of potential energy for
% palh2m2
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
% Datum: 2019-06-06 14:46
% Revision: 7254ec7b167830f9592b38d39d95d449e6fd98ef (2019-06-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh2m2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-06 14:45:19
% EndTime: 2019-06-06 14:45:19
% DurationCPUTime: 0.06s
% Computational Cost: add. (57->31), mult. (92->31), div. (0->0), fcn. (82->8), ass. (0->19)
t56 = sin(qJ(1));
t60 = cos(qJ(1));
t61 = g(1) * t60 + g(2) * t56;
t54 = sin(qJ(3));
t55 = sin(qJ(2));
t66 = t55 * pkin(4);
t67 = g(3) * (t54 * pkin(5) + t66);
t59 = cos(qJ(2));
t64 = t59 * pkin(4) + pkin(1);
t63 = g(3) * t66;
t62 = pkin(2) + t64;
t58 = cos(qJ(3));
t57 = cos(qJ(4));
t53 = sin(qJ(4));
t47 = g(1) * t56 - g(2) * t60;
t45 = t58 * pkin(5) + t62;
t44 = -t56 * t53 + t60 * t57;
t43 = t60 * t53 + t56 * t57;
t1 = [0, 0, 0, 0, 0, 0, -t61, t47, -g(3), 0, 0, 0, 0, 0, 0, 0, -g(3) * t55 - t61 * t59, -g(3) * t59 + t61 * t55, -t47, -t61 * pkin(1), 0, 0, 0, 0, 0, 0, -t61, -g(3), -t47, -t61 * t64 - t63, 0, 0, 0, 0, 0, 0, -g(3) * t54 - t61 * t58, -g(3) * t58 + t61 * t54, -t47, -t61 * t62 - t63, 0, 0, 0, 0, 0, 0, -t61, -g(3), -t47, -t61 * t45 - t67, 0, 0, 0, 0, 0, 0, -g(1) * t44 - g(2) * t43, g(1) * t43 - g(2) * t44, -g(3), -t67 - t61 * (pkin(3) + t45);];
U_reg  = t1;
