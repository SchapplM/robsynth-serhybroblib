% Calculate inertial parameters regressor of potential energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh2m1OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:27:44
% EndTime: 2020-05-03 00:27:44
% DurationCPUTime: 0.16s
% Computational Cost: add. (133->50), mult. (259->71), div. (0->0), fcn. (220->10), ass. (0->32)
t66 = sin(qJ(4));
t71 = cos(qJ(4));
t69 = sin(qJ(1));
t86 = g(2) * t69;
t74 = cos(qJ(1));
t87 = g(1) * t74;
t89 = -t87 - t86;
t54 = g(3) * t71 - t66 * t89;
t67 = sin(qJ(3));
t72 = cos(qJ(3));
t84 = t66 * g(3);
t79 = t71 * t89 + t84;
t51 = t54 * t72 - t67 * t79;
t68 = sin(qJ(2));
t73 = cos(qJ(2));
t90 = t67 * t54 + t79 * t72;
t49 = t51 * t68 + t90 * t73;
t88 = g(3) * pkin(5);
t85 = g(3) * t67;
t62 = pkin(6) * t66 + pkin(3);
t81 = pkin(4) * t86;
t78 = t72 * t89 + t85;
t77 = t51 * t73 - t68 * t90;
t76 = t89 * pkin(1) - t88;
t75 = pkin(4) * t84 - (pkin(4) * t71 + t62) * t87 - (g(3) * pkin(6) + t81) * t71 - t62 * t86;
t70 = cos(qJ(5));
t65 = sin(qJ(5));
t61 = t72 * pkin(3) + pkin(2);
t57 = -g(1) * t69 + g(2) * t74;
t55 = g(3) * t72 - t67 * t89;
t52 = -(-pkin(4) * t66 + pkin(6) * t71) * t87 + (g(3) * pkin(4) - pkin(6) * t86) * t71 + t66 * t81 + g(3) * t62;
t1 = [0, 0, 0, 0, 0, 0, t89, -t57, -g(3), -t88, 0, 0, 0, 0, 0, 0, g(3) * t68 + t89 * t73, g(3) * t73 - t68 * t89, -t57, t76, 0, 0, 0, 0, 0, 0, t55 * t68 + t78 * t73, t55 * t73 - t78 * t68, -t57, g(3) * (t68 * pkin(2) - pkin(5)) + t89 * (t73 * pkin(2) + pkin(1)), 0, 0, 0, 0, 0, 0, t49, t77, -t57, (pkin(3) * t85 + t61 * t89) * t73 + g(3) * (t61 * t68 - pkin(5)) - t89 * (t68 * t67 * pkin(3) - pkin(1)), 0, 0, 0, 0, 0, 0, t49 * t70 - t65 * t57, -t49 * t65 - t70 * t57, -t77, (pkin(2) * t89 + t52 * t67 + t75 * t72) * t73 + (g(3) * pkin(2) + t52 * t72 - t75 * t67) * t68 + t76;];
U_reg = t1;
