% Calculate inertial parameters regressor of potential energy for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh2m1DE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:29
% EndTime: 2020-05-02 23:52:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (61->30), mult. (109->41), div. (0->0), fcn. (90->8), ass. (0->22)
t65 = sin(qJ(1));
t76 = g(2) * t65;
t69 = cos(qJ(1));
t77 = g(1) * t69;
t79 = -t77 - t76;
t78 = pkin(5) * g(3);
t63 = sin(qJ(3));
t75 = t63 * g(3);
t67 = cos(qJ(3));
t60 = t67 * pkin(3) + pkin(2);
t64 = sin(qJ(2));
t74 = t64 * t60 - pkin(5);
t73 = t64 * t63 * pkin(3) - pkin(1);
t72 = -pkin(4) + t73;
t70 = t67 * t79 + t75;
t68 = cos(qJ(2));
t66 = cos(qJ(4));
t62 = sin(qJ(4));
t61 = pkin(3) * t75;
t55 = g(1) * t65 - g(2) * t69;
t54 = t67 * g(3) - t63 * t79;
t1 = [0, 0, 0, 0, 0, 0, t79, t55, -g(3), -t78, 0, 0, 0, 0, 0, 0, g(3) * t64 + t79 * t68, g(3) * t68 - t64 * t79, t55, t79 * pkin(1) - t78, 0, 0, 0, 0, 0, 0, t54 * t64 + t68 * t70, t54 * t68 - t64 * t70, t55, g(3) * (t64 * pkin(2) - pkin(5)) + t79 * (t68 * pkin(2) + pkin(1)), 0, 0, 0, 0, 0, 0, t79, g(3), t55, (t60 * t79 + t61) * t68 + t74 * g(3) - t79 * t73, 0, 0, 0, 0, 0, 0, t62 * t55 + t66 * t79, t55 * t66 - t62 * t79, -g(3), -(t60 * t68 - t72) * t77 + (-t60 * t76 + t61) * t68 + t72 * t76 + (-pkin(6) + t74) * g(3);];
U_reg = t1;
