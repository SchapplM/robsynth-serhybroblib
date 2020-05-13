% Implicit kinematic constraints of
% palh4m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:05
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = palh4m1IC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1IC_kinconstr_impl_mdh_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1IC_kinconstr_impl_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 23:05:39
% EndTime: 2020-04-11 23:05:40
% DurationCPUTime: 0.03s
% Computational Cost: add. (12->11), mult. (7->7), div. (0->0), fcn. (6->6), ass. (0->10)
unknown=NaN(3,1);
t1 = qJ(2) + qJ(4);
t2 = cos(t1);
t4 = cos(qJ(2));
t6 = sin(qJ(7));
t9 = sin(t1);
t11 = sin(qJ(2));
t13 = cos(qJ(7));
unknown(1,1) = t6 * pkin(1) - t2 * pkin(2) + t4 * qJ(3) - pkin(6) - pkin(7);
unknown(2,1) = -t13 * pkin(1) - t9 * pkin(2) + t11 * qJ(3);
unknown(3,1) = 0.5e1 / 0.2e1 * pi + qJ(2) + qJ(4) + qJ(8) - qJ(7);
h  = unknown;
