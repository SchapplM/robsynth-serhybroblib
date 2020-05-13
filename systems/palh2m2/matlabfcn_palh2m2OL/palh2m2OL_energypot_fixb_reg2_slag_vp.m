% Calculate inertial parameters regressor of potential energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh2m2OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 02:56:59
% EndTime: 2020-05-03 02:56:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (161->36), mult. (153->55), div. (0->0), fcn. (141->14), ass. (0->32)
t88 = cos(qJ(2));
t106 = -t88 * pkin(4) - pkin(1);
t85 = sin(qJ(1));
t89 = cos(qJ(1));
t105 = -g(1) * t89 - g(2) * t85;
t81 = qJ(2) + qJ(3);
t78 = qJ(4) + t81;
t77 = qJ(5) + t78;
t71 = sin(t77);
t102 = g(3) * t71;
t84 = sin(qJ(2));
t79 = t84 * pkin(4);
t82 = sin(qJ(6));
t101 = t85 * t82;
t86 = cos(qJ(6));
t100 = t85 * t86;
t99 = t89 * t82;
t98 = t89 * t86;
t97 = pkin(2) * sin(t81) + t79;
t75 = sin(t78);
t96 = pkin(5) * t75 + t97;
t95 = pkin(2) * cos(t81) - t106;
t83 = sin(qJ(3));
t87 = cos(qJ(3));
t91 = t83 * t88 + t84 * t87;
t90 = t83 * t84 - t87 * t88;
t76 = cos(t78);
t72 = cos(t77);
t69 = -g(1) * t85 + g(2) * t89;
t67 = pkin(5) * t76 + t95;
t66 = -g(3) * t72 - t105 * t71;
t1 = [0, 0, 0, 0, 0, 0, t105, -t69, -g(3), 0, 0, 0, 0, 0, 0, 0, -g(3) * t84 + t105 * t88, -g(3) * t88 - t105 * t84, t69, t105 * pkin(1), 0, 0, 0, 0, 0, 0, -t91 * g(3) - t105 * t90, t90 * g(3) - t105 * t91, t69, -g(3) * t79 - t105 * t106, 0, 0, 0, 0, 0, 0, -g(3) * t75 + t105 * t76, -g(3) * t76 - t105 * t75, t69, -g(3) * t97 + t105 * t95, 0, 0, 0, 0, 0, 0, t105 * t72 - t102, t66, t69, -g(3) * t96 + t105 * t67, 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t98 - t101) - g(2) * (t72 * t100 + t99) - t86 * t102, -g(1) * (-t72 * t99 - t100) - g(2) * (-t72 * t101 + t98) + t82 * t102, t66, -g(3) * (t71 * pkin(3) + t96) + t105 * (pkin(3) * t72 + t67);];
U_reg = t1;
